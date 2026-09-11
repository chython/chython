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
"""MDL reaction reading and writing.

The fixtures are the repository's own RDfiles.  `test/MR.rdf` earns its keep by being genuinely
mixed -- three V2000 CTABs and three V3000 ones in four records -- which is why the version sniff is
per-CTAB and not per-file.
"""

from chython.formats.ctfile import parse_v2000


_RXCTR_V2000 = ['reacting centre', '', '',
                '  2  1  0  0  0  0            999 V2000',
                '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                '    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
                # columns 16-18 are the bond topology, 19-21 the reacting centre status
                '  1  2  1  0  0  1  1',
                'M  END']


def test_reacting_centre_is_logged_not_silently_dropped():
    log = []
    ctab = parse_v2000(_RXCTR_V2000, log)
    assert ctab.bonds[0].reacting_center == 1
    _, _, log = ctab.build()
    assert any('reacting-centre' in x for x in log), log


def test_bond_topology_is_logged_as_unsupported():
    log = []
    ctab = parse_v2000(_RXCTR_V2000, log)
    assert ctab.bonds[0].topology == 1
    _, _, log = ctab.build()
    assert any(str(x).startswith('unsupported') and 'topology' in x for x in log), log


_RXN_V2000 = ['$RXN', 'ethanol to acetaldehyde', '', '',
              '  1  1',
              '$MOL', 'reactant', '', '',
              '  2  1  0  0  0  0            999 V2000',
              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
              '    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
              '  1  2  1  0  0  0  0',
              'M  END',
              '$MOL', 'product', '', '',
              '  2  1  0  0  0  0            999 V2000',
              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
              '    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
              '  1  2  2  0  0  0  0',
              'M  END']


def test_v2000_rxn_sides():
    from chython.formats.ctfile import parse_rxn

    log = []
    reaction = parse_rxn(_RXN_V2000, log)
    assert len(reaction.reactants) == 1
    assert len(reaction.products) == 1
    assert not reaction.agents
    assert reaction.title == 'ethanol to acetaldehyde'
    # the sides are not the same molecule: one C-O single, one C=O double
    assert next(iter(reaction.reactants[0].bonds())).order == 1
    assert next(iter(reaction.products[0].bonds())).order == 2


def test_v2000_third_count_is_read_and_logged_as_a_convention():
    """The counts line officially carries two fields.  A third is ubiquitous and unofficial."""
    from chython.formats.ctfile import parse_rxn

    lines = list(_RXN_V2000)
    lines[4] = '  1  1  1'
    lines += ['$MOL', 'agent', '', '',
              '  1  0  0  0  0  0            999 V2000',
              '    0.0000    0.0000    0.0000 Pd  0  0  0  0  0  0  0  0  0  0  0  0',
              'M  END']
    log = []
    reaction = parse_rxn(lines, log)
    assert len(reaction.agents) == 1
    assert any('third count' in x for x in log), log


def test_empty_rxn_is_an_empty_reaction_not_an_exception():
    """A `  0  0` counts line is an empty reaction, not an error: the record is not lost."""
    from chython.formats.ctfile import parse_rxn

    reaction = parse_rxn(['$RXN', '', '', '', '  0  0'], [])
    assert not reaction.reactants and not reaction.products and not reaction.agents


def test_missing_mol_block_logs_and_keeps_what_was_read():
    """The counts line promises two components and the file holds one."""
    from chython.formats.ctfile import parse_rxn

    lines = _RXN_V2000[:14]  # header, counts, and only the reactant block
    log = []
    reaction = parse_rxn(lines, log)
    assert len(reaction.reactants) == 1
    assert not reaction.products
    assert any('counts line promises' in x for x in log), log


def test_title_with_non_utf8_byte_round_trips_without_raising():
    """A non-UTF-8 name line arrives as lone surrogates from `surrogateescape` and stays that str: the
    reader stores it unchanged, so encoding it back with the same handler yields the original byte."""
    from chython.formats.ctfile import parse_rxn

    # U+DCE9 is the surrogateescape encoding of the byte 0xe9 (e.g. 'café' in latin-1).
    lines = ['$RXN', 'caf\udce9', '', '', '  0  0']
    reaction = parse_rxn(lines, [])
    assert reaction.title == 'caf\udce9'
    assert reaction.title.encode('utf8', 'surrogateescape') == b'caf\xe9'


def test_component_log_lines_are_attributed_and_not_deduplicated():
    """Two components triggering the same message produce two lines: a `x not in log` dedupe hides
    whether one or both had the defect, so merge_log prefixes `component N: `."""
    from chython.formats.ctfile import parse_rxn

    # A bond block that lists the same pair twice; Ctab.build logs 'duplicate bond 2 dropped'.
    _dup_mol = ['', '', '',
                '  2  2  0  0  0  0            999 V2000',
                '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                '    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
                '  1  2  1  0  0  0  0',
                '  1  2  2  0  0  0  0',  # duplicate pair -- will be dropped and logged
                'M  END']
    lines = ['$RXN', 'dup test', '', '', '  2  0',
             '$MOL'] + _dup_mol + ['$MOL'] + _dup_mol
    log = []
    reaction = parse_rxn(lines, log)
    assert len(reaction.reactants) == 2
    dup_lines = [x for x in log if 'duplicate bond' in x]
    assert len(dup_lines) == 2, f'expected 2 duplicate-bond log lines, got: {dup_lines}'
    assert any('component 1' in x for x in dup_lines), dup_lines
    assert any('component 2' in x for x in dup_lines), dup_lines


def test_v3000_rxn_from_the_repository_fixture(root):
    """MR.rdf record 4 is `$RXN V3000` with two inline CTABs and no `$MOL` delimiter."""
    from chython.formats.ctfile import V3000_STAMP, parse_rxn, sniff_rxn_version

    text = (root / 'test' / 'MR.rdf').read_text(encoding='utf8').split('\n')
    start = next(i for i, x in enumerate(text) if x.startswith('$RXN V3000'))
    end = next(i for i, x in enumerate(text[start:], start) if x.startswith('$DTYPE'))
    lines = [x.rstrip('\r') for x in text[start:end]]

    log = []
    assert sniff_rxn_version(lines, log) == V3000_STAMP
    reaction = parse_rxn(lines, log)
    assert len(reaction.reactants) == 1, log
    assert len(reaction.products) == 1, log
    assert reaction.reactants[0].atom_count == 6
    assert reaction.title == 'title6'


def test_v3000_rxn_agents_block():
    from chython.formats.ctfile import parse_rxn
    from chython.core._core import element_symbols

    lines = ['$RXN V3000', 'with an agent', '', '',
             'M  V30 COUNTS 1 1 1',
             'M  V30 BEGIN REACTANT',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 1 0 0 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 C 0 0 0 0',
             'M  V30 END ATOM',
             'M  V30 END CTAB',
             'M  V30 END REACTANT',
             'M  V30 BEGIN PRODUCT',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 1 0 0 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 O 0 0 0 0',
             'M  V30 END ATOM',
             'M  V30 END CTAB',
             'M  V30 END PRODUCT',
             'M  V30 BEGIN AGENT',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 1 0 0 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 Pd 0 0 0 0',
             'M  V30 END ATOM',
             'M  V30 END CTAB',
             'M  V30 END AGENT',
             'M  END']
    log = []
    reaction = parse_rxn(lines, log)
    assert [element_symbols()[m.atom(next(iter(m.atom_numbers))).element]
            for m in reaction.agents] == ['Pd'], log
    # V3000 spells agents, so nothing about a third count is logged here
    assert not any('third count' in x for x in log), log


def test_v3000_truncated_ctab_logs_cause_and_consequence():
    """A CTAB missing END CTAB: parse_ctab_block logs the cause and the skip walk the consequence.
    Two layers, two claims, so silencing one must not silence the other."""
    from chython.formats.ctfile import parse_rxn

    lines = ['$RXN V3000', 'trunc test', '', '',
             'M  V30 BEGIN REACTANT',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 1 0 0 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 C 0 0 0 0',
             'M  V30 END ATOM',
             # END CTAB omitted intentionally
             'M  V30 END REACTANT',
             'M  V30 BEGIN PRODUCT',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 1 0 0 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 O 0 0 0 0',
             'M  V30 END ATOM',
             'M  V30 END CTAB',
             'M  V30 END PRODUCT',
             'M  END']
    log = []
    reaction = parse_rxn(lines, log)
    assert len(reaction.reactants) == 1, log
    # The product CTAB is unreachable once the truncated CTAB swallows the rest of the body.
    assert len(reaction.products) == 0, log
    # The cause: parse_ctab_block found no END CTAB.
    assert any('END CTAB missing' in x and 'end of input' in x for x in log), log
    # The consequence: the skip walk confirms components after this point are lost.
    assert any('END CTAB missing' in x and 'after this point' in x for x in log), log


def test_emit_rxn_v2000_shape():
    from chython.formats.ctfile import V2000_STAMP, emit_rxn, parse_rxn

    reaction = parse_rxn(_RXN_V2000, [])
    lines, log = emit_rxn(reaction, version=V2000_STAMP, title='rewritten')
    assert lines[0] == '$RXN'
    assert lines[1] == 'rewritten'
    assert lines[4] == '  1  1'
    assert lines.count('$MOL') == 2
    assert lines[-1] == 'M  END'


def test_emit_rxn_v2000_agents_use_the_third_count_and_log_it():
    """Agents on a V2000 write do NOT escalate to V3000 -- the third count is what tools read."""
    from chython.formats.ctfile import V2000_STAMP, emit_rxn, parse_rxn

    lines = list(_RXN_V2000)
    lines[4] = '  1  1  1'
    lines += ['$MOL', 'agent', '', '',
              '  1  0  0  0  0  0            999 V2000',
              '    0.0000    0.0000    0.0000 Pd  0  0  0  0  0  0  0  0  0  0  0  0',
              'M  END']
    reaction = parse_rxn(lines, [])
    out, log = emit_rxn(reaction, version=V2000_STAMP)
    assert out[4] == '  1  1  1'
    assert out.count('$MOL') == 3
    assert any('third count' in x for x in log), log


def test_emit_rxn_v3000_nests_bare_ctabs():
    from chython.formats.ctfile import V3000_STAMP, emit_rxn, parse_rxn

    reaction = parse_rxn(_RXN_V2000, [])
    lines, log = emit_rxn(reaction, version=V3000_STAMP)
    assert lines[0] == '$RXN V3000'
    assert lines[4] == 'M  V30 COUNTS 1 1'
    assert 'M  V30 BEGIN REACTANT' in lines
    assert 'M  V30 END PRODUCT' in lines
    assert '$MOL' not in lines
    # 'reactant' is the first component's molfile title: the component headers must not survive
    assert 'reactant' not in lines
    # exactly one M  END, at the end of the record
    assert lines.count('M  END') == 1 and lines[-1] == 'M  END'


def test_rxn_round_trips_both_versions():
    from chython.formats.ctfile import V2000_STAMP, V3000_STAMP, emit_rxn, parse_rxn

    original = parse_rxn(_RXN_V2000, [])
    for version in (V2000_STAMP, V3000_STAMP):
        lines, _ = emit_rxn(original, version=version)
        again = parse_rxn(lines, [])
        assert len(again.reactants) == len(original.reactants)
        assert len(again.products) == len(original.products)
        assert ([m.atom_count for m in again.molecules()]
                == [m.atom_count for m in original.molecules()])


def test_emit_rxn_writes_the_name_line_byte_for_byte():
    """The RXN writer takes NO loss on a name line, so there is nothing left for it to report.

    Two tests went with the loss: one on the replacement count in the `unsupported: ` line, one on a
    name whose real bytes encode U+FFFD staying quiet.  Both had the removed strict-decode detector as
    their whole subject, and `chython/formats/test/test_log_prefix.py` guards the prefix constant for
    every message that remains.
    """
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import V2000_STAMP, V3000_STAMP, emit_rxn

    rxn = ReactionContainer(title=b'caf\xe9 reaction')   # \xe9 is not UTF-8 on its own
    for version in (V2000_STAMP, V3000_STAMP):
        lines, log = emit_rxn(rxn, version=version)
        assert not log, log
        assert lines[1].encode('utf8', 'surrogateescape') == b'caf\xe9 reaction'

    # A caller-supplied title overrides the stored one and is written as given.
    lines, log = emit_rxn(rxn, version=V2000_STAMP, title='plain title')
    assert not log and lines[1] == 'plain title'


def test_emit_rxn_keeps_the_byte_and_the_stream_is_what_encodes_it(tmp_path):
    """INVERTED: this asserted the output carried no lone surrogate.  Keeping the byte is the point now.

    A file opened WITHOUT `errors=` is strict UTF-8 and refuses the line; `_FileBacked._open` opens with
    `surrogateescape`, which is why the library's own writers put the original byte back.
    """
    from pytest import raises
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import V2000_STAMP, emit_rxn

    rxn = ReactionContainer(title=b'caf\xe9 reaction')  # \xe9 is not valid UTF-8 alone
    lines, log = emit_rxn(rxn, version=V2000_STAMP)
    assert not log and lines[1] == 'caf\udce9 reaction'
    # `encoding` and `newline` stated on both writes, because the assertion is about bytes: the default
    # codec is the locale's and the default newline is the host's, so the same call writes `caf\xe9` under
    # one and `caf?` or a CRLF under another.  This is the pair `_FileBacked._open` uses.
    with raises(UnicodeEncodeError):
        (tmp_path / 'strict.rxn').write_text('\n'.join(lines), encoding='utf-8')
    (tmp_path / 'ok.rxn').write_text('\n'.join(lines), encoding='utf-8', errors='surrogateescape',
                                     newline='\n')
    assert (tmp_path / 'ok.rxn').read_bytes().split(b'\n')[1] == b'caf\xe9 reaction'


def test_emit_rxn_v2000_over_999_components_raises():
    """A V2000 counts field is three characters wide, so 1000 components shift every later field:
    MalformedCtfile, naming V3000 as the fix."""
    from chython.core._core import MoleculeContainer
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import V2000_STAMP, emit_rxn
    from chython.formats.ctfile._errors import MalformedCtfile

    m = MoleculeContainer()
    with m.edit() as e:
        e.add_atom(6)
    rxn = ReactionContainer([m] * 1000, [m], [])
    try:
        emit_rxn(rxn, version=V2000_STAMP)
        assert False, 'expected MalformedCtfile'
    except MalformedCtfile as exc:
        assert 'V3000' in str(exc), str(exc)


def test_v3000_nested_role_opening_is_logged():
    """V3000 requires explicit END REACTANT/PRODUCT/AGENT framing, so a BEGIN PRODUCT before the
    preceding END is reported -- the following CTABs are still attributed to the new role."""
    from chython.formats.ctfile import parse_rxn

    lines = ['$RXN V3000', 'nested role', '', '',
             'M  V30 BEGIN REACTANT',
             # END REACTANT omitted -- next line opens PRODUCT while REACTANT is still open
             'M  V30 BEGIN PRODUCT',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 1 0 0 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 O 0 0 0 0',
             'M  V30 END ATOM',
             'M  V30 END CTAB',
             'M  V30 END PRODUCT',
             'M  END']
    log = []
    reaction = parse_rxn(lines, log)
    assert len(reaction.products) == 1, log
    assert any('role framing is broken' in x for x in log), log


def test_parse_rxn_record_returns_a_reaction_carrying_its_meta():
    """An RDfile ``$RFMT`` record's ``$DTYPE`` pairs land on `ReactionContainer.meta`; the framing
    facts no container holds come out of ``header=``.  The name line is not among them -- that is
    ``reaction.title``.
    """
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import parse_rxn_record

    header = {}
    reaction = parse_rxn_record(_RXN_V2000, {'TEMP': '100'}, header=header)
    assert isinstance(reaction, ReactionContainer) and reaction.meta == {'TEMP': '100'}
    assert reaction.title == 'ethanol to acetaldehyde'
    assert set(header) == {'version', 'program', 'comment'}
