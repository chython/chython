# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""What survives a round trip when the caller keeps **the molecule alone**.

The title, the S-groups and the aliases live on the molecule, not beside it as reader state, so
``for mol in SDFRead(f): w.write(mol)`` preserves them.  Fidelity is the default: the writers still
take ``title=`` and ``sgroups=`` and those still override, but passing nothing means "what the
molecule says" rather than "empty".
"""
from contextlib import contextmanager
from io import StringIO
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp

from pytest import fixture, mark

# The one definition of "run something under chython 2", shared with the core's differentials.
from chython.core.test import oracle

from .._sdf import emit_record, parse_record
from .._sgroup import NO_INDEX, SGroupStore
from .._facade import mol
from .._stream import SDFRead, SDFWrite
from .._v3000 import V3000_STAMP


# Phenol, carrying one of everything the fidelity promise covers: a title, a DAT S-group with a
# FIELDDISP anchor, an SRU whose keyword is not interpreted, an atom alias, two SDF data fields.
RICH = '''phenol
  probe
comment line
  7  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.2124    0.7000    0.0000 C   0  0
    2.4249    0.0000    0.0000 C   0  0
    2.4249   -1.4000    0.0000 C   0  0
    1.2124   -2.1000    0.0000 C   0  0
    0.0000   -1.4000    0.0000 C   0  0
   -1.2124    0.7000    0.0000 O   0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
A    7
OH-label
M  STY  2   1 DAT   2 SRU
M  SAL   1  1   7
M  SDT   1 pKa
M  SDD   1     1.2000    3.4000    DAU   ALL  1       6
M  SED   1 9.95
M  SAL   2  2   3   4
M  SMT   2 n
M  END
>  <SOURCE>
public reference compound

>  <NOTE>
a second field

$$$$
'''


@fixture
def header():
    """`parse_record`'s ``header=`` out-parameter, filled by the `molecule` fixture below.

    A fixture of its own so a test can ask for the store without re-reading the record: pytest hands
    both fixtures the same dict within one test.
    """
    return {}


@fixture
def molecule(header):
    return parse_record(RICH.split('\n'), header=header)


# the molecule carries it

def test_the_molecule_carries_the_title(molecule):
    """`mol.title`, not just `reader.title`.  A title on the reader is lost by `[m for m in f]`."""
    assert molecule.title == 'phenol'


def test_a_non_utf8_name_line_round_trips_byte_for_byte(tmp_path):
    """THE WHOLE JUSTIFICATION FOR `title` BEING `str`, end to end through a real file.

    `\xe9` alone is not valid UTF-8.  The reader opens the file, so the `surrogateescape` decode is the
    library's own, and the writer puts the same byte back.  The substituted name line carries no trailing
    whitespace: `_v2000.py` rstrips it, which is the format's rule and not this promise's business.
    """
    source = tmp_path / 'in.sdf'
    source.write_bytes(RICH.replace('phenol', 'NAME', 1).encode('utf8')
                       .replace(b'NAME', b'caf\xe9', 1))
    with SDFRead(source) as f:
        mol = next(iter(f))
    assert isinstance(mol.title, str) and mol.title == 'caf\udce9'

    out = tmp_path / 'out.sdf'
    with SDFWrite(out) as w:
        w.write(mol)
    assert out.read_bytes().split(b'\n', 1)[0] == b'caf\xe9', 'the name line is the byte the file held'


def test_the_molecule_carries_the_sgroups(molecule):
    mol = molecule
    assert len(mol.sgroups) == 2
    assert [s['type'] for s in mol.sgroups] == [b'DAT', b'SRU']


def test_the_molecule_carries_the_alias(molecule):
    """`A  <n>` lines.  Stored as single-atom S-group records, exposed as {n: bytes}."""
    mol = molecule
    assert list(mol.aliases.values()) == [b'OH-label']
    # An alias is not reported as an S-group record, or a consumer counting DAT groups is off by
    # the number of aliases in the file.
    assert len(mol.sgroups) == 2


def test_the_dat_group_keeps_its_name_data_and_display_anchor(molecule):
    dat, = [s for s in molecule.sgroups if s['type'] == b'DAT']
    assert dat['name'] == b'pKa'
    assert dat['data'] == (b'9.95',)
    assert dat['disp'] == (1.2, 3.4)
    assert b'DAU' in dat['disp_tail']


def test_an_uninterpreted_sgroup_keyword_survives_on_the_molecule(molecule):
    """`SRU`'s `SMT` is not modelled and must still be on the molecule: what survives is not
    limited to what this library understands.
    """
    sru, = [s for s in molecule.sgroups if s['type'] == b'SRU']
    assert (b'LABEL', b'n') in sru['fields']


def test_sgroup_atom_references_are_stable_ids_not_positions(molecule):
    """A reference that survives a remap is the only kind worth storing."""
    mol = molecule
    dat, = [s for s in mol.sgroups if s['type'] == b'DAT']
    # `Atom.element` is the atomic number as an `int`, not the symbol -- so oxygen is 8.
    oxygen, = [n for n in mol.atom_numbers if mol.atom(n).element == 8]
    assert dat['atoms'] == (oxygen,)


# the store <-> molecule bridge

def test_store_round_trips_through_the_molecule(molecule, header):
    """`to_molecule` then `from_molecule` is the identity on everything the record models."""
    before = header['sgroups']
    after = SGroupStore.from_molecule(molecule)

    assert [s.type for s in after.records] == [s.type for s in before.records]
    assert [s.name for s in after.records] == [s.name for s in before.records]
    assert [s.data for s in after.records] == [s.data for s in before.records]
    assert [s.atoms for s in after.records] == [s.atoms for s in before.records]
    assert [s.bonds for s in after.records] == [s.bonds for s in before.records]
    assert [s.disp for s in after.records] == [s.disp for s in before.records]
    assert [s.fields for s in after.records] == [s.fields for s in before.records]
    assert after.aliases == before.aliases


def test_sgroup_numbers_round_trip_verbatim(molecule):
    """`index` is the file's number, not ours, and `NO_INDEX` must not collapse into 0.

    The sentinel is 0xFFFF because 0 is a number a file may legitimately carry, so "unnumbered" and
    "numbered 0" stay distinguishable across the storage.
    """
    after = SGroupStore.from_molecule(molecule)
    assert [s.index for s in after.records] == [1, 2]
    assert all(s.parent == NO_INDEX for s in after.records)


# fidelity by default

def test_writing_a_molecule_with_no_arguments_keeps_its_title_and_sgroups(molecule):
    """`write(mol)` -- nothing else -- must not produce an untitled record with no S-groups."""
    buffer = StringIO()
    with SDFWrite(buffer) as f:
        f.write(molecule)
    text = buffer.getvalue()

    assert text.split('\n')[0] == 'phenol'
    assert 'M  STY' in text
    assert 'pKa' in text and '9.95' in text


def test_the_naive_read_write_loop_preserves_everything():
    """Read, write, read again, keeping only the molecule between the two."""
    buffer = StringIO()
    with SDFRead(StringIO(RICH)) as r, SDFWrite(buffer) as w:
        for mol in r:
            w.write(mol)

    again = parse_record(buffer.getvalue().split('\n'))
    assert again.title == 'phenol'
    assert again.meta == {'SOURCE': 'public reference compound', 'NOTE': 'a second field'}
    assert [s['type'] for s in again.sgroups] == [b'DAT', b'SRU']
    assert list(again.aliases.values()) == [b'OH-label']
    dat, = [s for s in again.sgroups if s['type'] == b'DAT']
    assert dat['name'] == b'pKa' and dat['data'] == (b'9.95',)


@mark.parametrize('version', [None, V3000_STAMP])
def test_both_emitters_default_to_the_molecules_own(molecule, version):
    """V3000 must not be the version where fidelity silently stops."""
    kwargs = {} if version is None else {'version': version}
    lines, _ = emit_record(molecule, **kwargs)
    text = '\n'.join(lines)
    assert 'phenol' in text
    assert 'pKa' in text


def test_an_explicit_title_still_overrides(molecule):
    """Fidelity is the default, not a lock: re-titling a record stays possible."""
    lines, _ = emit_record(molecule, title='renamed')
    assert lines[0] == 'renamed'


def test_an_explicit_empty_sgroup_store_still_suppresses(molecule):
    """`sgroups=SGroupStore()` means "write none" and must not fall back to the molecule's: "I did
    not say" and "I said none" are different, and conflating them makes suppression impossible.
    """
    lines, _ = emit_record(molecule, sgroups=SGroupStore())
    assert 'M  STY' not in '\n'.join(lines)


# a writer does not touch its input

def test_writing_does_not_mutate_the_molecule(molecule):
    """IO is not a mutator of representation, in the writing direction too: a serialiser that
    repairs, re-lays-out or kekulises its argument changes the caller's data.  The canonical bytes are
    the cheapest total statement of "nothing changed".
    """
    mol = molecule
    before = (mol.title, mol.sgroups, mol.aliases, bytes(mol.canonical_bytes),
              [mol.xy_of(n) for n in mol.atom_numbers])

    emit_record(mol)
    emit_record(mol, version=V3000_STAMP)

    assert (mol.title, mol.sgroups, mol.aliases, bytes(mol.canonical_bytes),
            [mol.xy_of(n) for n in mol.atom_numbers]) == before


def test_reading_does_not_normalise(molecule):
    """An aromatic input stays aromatic; a Kekule input stays Kekule.

    `RICH` draws phenol Kekule, so the molecule comes back with zero aromatic bonds; a reader running
    `thiele()` for the caller would fail here.
    """
    assert molecule.aromatic_bond_count == 0

    aromatic = RICH.replace('  1  2  2  0', '  1  2  4  0')
    assert mol(aromatic).aromatic_bond_count == 1


# `clean_stereo()` survives the trip through a molfile

# L-alanine, with an up wedge from the alpha carbon to its nitrogen.  The configuration is stated
# only by the drawing, so the wedge is all the reader has to work from.
WEDGED = '''L-alanine
  probe
wedge on the alpha carbon
  6  5  0  0  1  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
   -1.2124    0.7000    0.0000 C   0  0
    0.0000    1.4000    0.0000 N   0  0
    1.2124   -0.7000    0.0000 C   0  0
    1.2124   -2.1000    0.0000 O   0  0
    2.4249    0.0000    0.0000 O   0  0
  1  2  1  0
  1  3  1  1
  1  4  1  0
  4  5  2  0
  4  6  1  0
M  END
$$$$
'''


def _configured(mol):
    return {n: mol.parity_of(n) for n in mol.atom_numbers if mol.parity_of(n)}


def test_the_wedge_really_is_the_only_statement_of_the_configuration():
    """The control: reading WEDGED configures the alpha carbon, so the round trip below is not
    running on a record that never carried a configuration.  The reader derives the parity from the
    wedge, which is what makes the next test non-trivial.
    """
    m = mol(WEDGED)
    assert len(_configured(m)) == 1, 'the alpha carbon, and nothing else'
    assert m.wedges(), 'and the wedge it was derived from is stored'

    lines, _ = emit_record(m)
    assert _configured(mol('\n'.join(lines))) == _configured(m), \
        'an untouched molecule keeps its configuration across a write and a read'


def test_a_cleaned_molecule_does_not_get_its_parities_back_from_its_own_wedges():
    """The wedge half of `clean_stereo`, which only a round trip can pin.

    `wedges_for_write` returns the wedges a molecule already carries rather than re-deriving them,
    and `assign_parities` reads parities back out of a drawing.  So a `clean_stereo` that cleared the
    parities and left the wedges is invisible in the arena and undone by one molfile.
    """
    m = mol(WEDGED)
    report = m.clean_stereo()
    assert 'parities' in report and 'wedges' in report, 'both kinds were there to be cleared'
    assert _configured(m) == {} and m.wedges() == []

    lines, log = emit_record(m)
    again = mol('\n'.join(lines))
    assert again.wedges() == [], 'nothing was drawn, so there was nothing to read a sign out of'
    assert _configured(again) == {}, 'and no parity came back'
    # otherwise intact: same constitution, same layout
    assert bytes(again.canonical_bytes) == bytes(m.canonical_bytes)
    assert [again.xy_of(n) for n in again.atom_numbers] == [m.xy_of(n) for n in m.atom_numbers]


# what this library does not understand, both ways

UNKNOWN_KEYWORD = '''keyword probe


  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 O 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 1 ATOMS=(1 2) LABEL=Et NATREPLACE=SOME/THING
M  V30 END SGROUP
M  V30 END CTAB
M  END
$$$$
'''


def test_an_unknown_v3000_sgroup_keyword_round_trips_verbatim():
    """A keyword nothing here interprets survives: `NATREPLACE` has no slot, no validation and no
    meaning attached, and must still come back out of the writer.  A refactor that starts filtering
    `fields` to known keywords fails here.
    """
    m = mol(UNKNOWN_KEYWORD)
    sup, = m.sgroups
    assert (b'NATREPLACE', b'SOME/THING') in sup['fields']

    lines, log = emit_record(m, version=V3000_STAMP)
    assert 'NATREPLACE=SOME/THING' in '\n'.join(lines)
    assert not log, log


CSTATE_VECTOR = UNKNOWN_KEYWORD.replace('LABEL=Et NATREPLACE=SOME/THING',
                                        'XBONDS=(1 1) LABEL=Et CSTATE=(4 1 0.808958 -0.158697 0)')


def test_a_cstate_vector_crosses_the_arena_byte_for_byte():
    """The attachment vector is a string in the model and bytes in the arena, and `_to_dict` is where
    it is encoded.

    Every other CSTATE test hands the writer the READER'S store, so the vector never crosses into the
    molecule; this one goes through `mol()`, which is the path a caller takes.  What the encoding
    prevents: the blob run copies `len(item)` bytes from a `bytes` payload pointer, so an unencoded
    `str` stored sixteen bytes of CPython object header and a truncated vector, and the writer then
    raised `TypeError` on the bytes it read back.  `structure_put_blob` now refuses a non-`bytes` item
    outright, so the encoding cannot be forgotten again in silence.
    """
    m = mol(CSTATE_VECTOR)
    sup, = m.sgroups
    assert sup['cstates'] == (((1, 2), b'0.808958 -0.158697 0'),)

    text = mol(m, version=3000)
    assert 'CSTATE=(4 1 0.808958 -0.158697 0)' in text.replace(' -\nM  V30 ', ' ')
    again, = mol(text).sgroups
    assert again['cstates'] == sup['cstates'], 'the vector moved on the second crossing'


HYPHEN_LABEL = UNKNOWN_KEYWORD.replace('LABEL=Et NATREPLACE=SOME/THING', 'LABEL="NH3+Cl-"') \
    .replace('M  V30 END SGROUP', 'M  V30 END SGROUP\nM  V30 BEGIN COLLECTION'
             '\nM  V30 MDLV30/STEABS ATOMS=(1 1)\nM  V30 END COLLECTION')


def test_a_label_ending_in_a_hyphen_does_not_swallow_the_rest_of_the_ctab():
    """The blocks after the S-group must still be there, which is what makes the quoting load-bearing
    rather than cosmetic.

    Bare, the label's own hyphen ends the physical line and the reader joins `END SGROUP` onto it --
    so the S-group block never closes, the collection block after it is read as S-group content, and
    the enhanced stereo is gone with no diagnostic naming the label.  A common salt drawing:
    ChemDraw writes ``LABEL="NH3+Cl-"`` and the writer has to as well.
    """
    m = mol(HYPHEN_LABEL)
    sup, = m.sgroups
    assert sup['fields'] == ((b'LABEL', b'NH3+Cl-'),)
    assert m.stereo_groups() == {(1, 0): [1]}

    log = []
    again = mol(mol(m, version=3000), log=log)
    assert not log, log
    assert again.stereo_groups() == m.stereo_groups(), 'the collection was read as S-group content'
    assert again.sgroups[0]['fields'] == sup['fields']


def test_an_unrecognised_v2000_property_line_is_dropped_and_says_so():
    """A declared gap: an unrecognised V2000 property line is not preserved.

    The molecule has segments for a title, S-groups and aliases and none for an opaque line, so
    closing the gap needs a core segment rather than a parser change.  The loss is reported per
    record, which a caller sweeping 40,000 records can act on.  When the segment exists this test
    fails on its first assertion, which is how it becomes obsolete.
    """
    text = RICH.replace('M  END', 'M  ZZZ  1 SOMETHING\nM  END')
    molecule = parse_record(text.split('\n'))

    assert any('M  ZZZ' in str(x) and 'dropped' in str(x) for x in molecule.log), molecule.log
    assert 'ZZZ' not in '\n'.join(emit_record(molecule)[0])
    # The gap is confined to the unknown line: a record with one junk property still has its
    # S-groups.
    assert molecule.title == 'phenol'
    assert len(molecule.sgroups) == 2


# the differential oracle: chython 2.24

def _oracle(source):
    """Run `source` under the pinned, isolated chython 2.24 and hand back its stdout.

    A subprocess against an installed 2.24, not an import: V2 is not in this tree.  The spawn lives
    in `chython/core/test/oracle.py`, the only place that starts an oracle interpreter; it passes
    `-I`, checks the version and checks which `chython` the child imported, on every call.
    """
    return oracle.ask_text(source)


@contextmanager
def _record_file():
    """`RICH` on disk, for the oracle to open.  Passed by path so the record text stays in one place."""
    path = Path(mkdtemp()) / 'probe.mol'
    path.write_text(RICH)
    try:
        yield path
    finally:
        rmtree(path.parent, ignore_errors=True)


def test_the_oracle_is_the_pinned_version_and_is_not_this_tree():
    """Pinned, because "whatever chython 2 is installed" is not a fixed comparison, and separate,
    because an oracle importing `./chython/` would agree with this parser about everything.  Both
    checks run on every `_oracle` call; their own tests are in `chython/core/test/test_oracle.py`.
    """
    oracle.require()
    oracle.verify()
    assert oracle.probe()[0] == oracle.VERSION


def test_v3_reads_the_same_constitution_as_chython_2():
    """The parser must agree with V2 about the constitution.

    Not compared by SMILES string: the two versions' canonical writers are independent and are
    allowed to disagree about output.
    """
    with _record_file() as path:
        out = _oracle('''
from chython import mdl_mol
m = mdl_mol(open(%r).read())
print('ATOMS', m.atoms_count)
print('BONDS', m.bonds_count)
print('SMILES', str(m))
''' % str(path))

    ours = mol(RICH)
    assert 'ATOMS %d' % ours.atom_count in out
    assert 'BONDS %d' % ours.bond_count in out
    # The two canonical SMILES writers are independent, so the strings may differ; the molecule may
    # not -- same heavy-atom composition, same aromatic/Kekule representation as drawn.
    v2 = out.split('SMILES', 1)[1].strip()
    assert v2.count('O') == 1 and v2.count('C') == 6
    assert ours.aromatic_bond_count == 0 and v2.islower() is False


def test_chython_2_loses_every_sgroup_and_alias():
    """A deliberate divergence, pinned: chython 2.24 has no S-group storage, logs `ignored line` /
    `ignored data` for the alias and the DAT group, and returns a molecule with neither a `sgroups`
    nor an `aliases` attribute.  V2 is the oracle for the constitution only.  If a future chython 2
    gained S-group support this fails, and the response is to compare rather than to delete it.
    """
    with _record_file() as path:
        out = _oracle('''
from chython import mdl_mol
m = mdl_mol(open(%r).read())
print('SGROUPS', hasattr(m, 'sgroups'))
print('ALIASES', hasattr(m, 'aliases'))
print('LOG', 'ignored' in repr(dict(m.meta)))
''' % str(path))

    assert 'SGROUPS False' in out, 'chython 2 grew S-group storage; compare rather than assume'
    assert 'ALIASES False' in out
    assert 'LOG True' in out, 'chython 2 no longer even logs what it drops'

    # the divergence, asserted on our side
    ours = mol(RICH)
    assert len(ours.sgroups) == 2
    assert list(ours.aliases.values()) == [b'OH-label']
