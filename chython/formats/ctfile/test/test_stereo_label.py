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
"""``STEREOLABEL``: a configuration word in the field the spec reserves for arbitrary data.

A large body of records states the partition this way and no other, V2000 having no syntax for an
enhanced-stereo collection at all.  Storing that word where no query can reach it is not reading it, so
the reader promotes it -- under two rules, both of which these tests state:

* a collection the record DOES state wins, the spec-defined statement being the one to believe;
* a promotion never invents a parity, so a label on an atom no configuration reaches is logged and
  dropped rather than turned into a member with nothing inside it.

The fixtures are two public compounds written by chython's own V2000 writer, because what is under test
is the ``DAT`` record and the promotion, not a vendor's idea of a column.
"""

from pytest import mark

from .._sgroup import STEREOLABEL, promote_stereo_labels, stereo_labels
from .._v2000 import emit_v2000, parse_v2000
from .._v3000 import parse_v3000
from ....core import INFO, LOST, REPAIRED, STEREO_ABS, STEREO_AND


#: L-alanine, one tetrahedral centre at atom 2, drawn with a down wedge to atom 1.
_ALANINE = ['L-alanine', '  chython', '',
            '  6  5  0  0  1  0            999 V2000',
            '    0.3572    1.0312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
            '    0.3572    0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
            '    1.0717   -0.2062    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0',
            '   -0.3572   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
            '   -0.3572   -1.0312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
            '   -1.0717    0.2063    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
            '  2  1  1  6  0  0  0',
            '  2  3  1  0  0  0  0',
            '  2  4  1  0  0  0  0',
            '  4  5  2  0  0  0  0',
            '  4  6  1  0  0  0  0']

#: (2R,3R)-butane-2,3-diol: TWO centres, at atoms 2 and 4, so one label can cover a pair and two
#: labels can ask for two collections.
_DIOL = ['(2R,3R)-butane-2,3-diol', '  chython', '',
         '  6  5  0  0  1  0            999 V2000',
         '    0.3572    1.0312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
         '    0.3572    0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
         '    1.0717   -0.2062    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
         '   -0.3572   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
         '   -0.3572   -1.0312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
         '   -1.0717    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
         '  2  1  1  6  0  0  0',
         '  2  3  1  0  0  0  0',
         '  2  4  1  0  0  0  0',
         '  4  5  1  0  0  0  0',
         '  4  6  1  6  0  0  0']

#: The same alanine as a V3000 record, which is the only version that can also state a collection --
#: `STERAC1` on the centre, so the two statements are side by side and one of them has to win.
_ALANINE_V3000 = ['L-alanine', '  chython', '',
                  '  0  0  0  0  0  0            999 V3000',
                  'M  V30 BEGIN CTAB',
                  'M  V30 COUNTS 6 5 0 0 0',
                  'M  V30 BEGIN ATOM',
                  'M  V30 1 C 0.3572 1.0312 0 0',
                  'M  V30 2 C 0.3572 0.2062 0 0',
                  'M  V30 3 N 1.0717 -0.2062 0 0',
                  'M  V30 4 C -0.3572 -0.2062 0 0',
                  'M  V30 5 O -0.3572 -1.0312 0 0',
                  'M  V30 6 O -1.0717 0.2063 0 0',
                  'M  V30 END ATOM',
                  'M  V30 BEGIN BOND',
                  'M  V30 1 1 2 1 CFG=3',
                  'M  V30 2 1 2 3',
                  'M  V30 3 1 2 4',
                  'M  V30 4 2 4 5',
                  'M  V30 5 1 4 6',
                  'M  V30 END BOND']


def _dat(no, atoms, data, name=STEREOLABEL):
    """One ``DAT`` group's four V2000 property lines: type, atom list, field name, datum."""
    return ['M  STY  1 %3d DAT' % no,
            'M  SAL %3d%3d' % (no, len(atoms)) + ''.join('%4d' % a for a in atoms),
            'M  SDT %3d ' % no + f'{name:<30s}',
            'M  SED %3d %s' % (no, data)]


def _read(*groups, lines=_ALANINE, ignore_stereo=False):
    """`lines` plus one ``DAT`` group per `(atoms, data)` or `(atoms, data, name)` triple."""
    properties = []
    for no, group in enumerate(groups, 1):
        properties.extend(_dat(no, *group[:1], *group[1:]))
    ctab = parse_v2000(lines + properties + ['M  END'], [])
    return ctab.build(ignore_stereo=ignore_stereo)


def _rules(molecule):
    return [x.rule for x in molecule.log]


# ------------------------------------------------------------------ stereo_labels(): the accessor

@mark.parametrize('datum,label,relative', [
    ('R', 'R', False),
    ('(R)', 'R', False),               # the parentheses are punctuation
    ('*R', 'R', True),                 # the marker is lifted out, not left in the string
    ('(*R)', 'R', True),
    ('*(R)', 'R', True),
    ('R*', 'R', True),                 # the same word with the marker on the other side
    ('(R*)', 'R', True),
    ('(R)*', 'R', True),
    ('*S*', 'S', True),
    ('  (S)  ', 'S', False),
    ('r', 'R', False),                 # one word, whatever case it arrived in
    ('rs', 'RS', False),
    ('*RS', 'RS', True),
    ('E', 'E', False),                 # read, and not promoted -- see below
    ('&', '&', False),
])
def test_the_payload_is_normalised_three_ways_and_the_rest_is_left_alone(datum, label, relative):
    molecule = _read(((2,), datum))[0]
    assert stereo_labels(molecule) == [((2,), (), label, relative)]


def test_the_field_name_matches_without_case():
    """``FIELDNAME`` is a key a producer chose, not a spelling the spec fixes, and one body of records
    writes both spellings."""
    for name in ('STEREOLABEL', 'StereoLabel', 'stereolabel'):
        molecule = _read(((2,), '(R)', name))[0]
        assert stereo_labels(molecule) == [((2,), (), 'R', False)], name


def test_another_field_name_is_not_a_stereo_label():
    molecule = _read(((2,), '(R)', 'COMMENT'))[0]
    assert stereo_labels(molecule) == []
    assert molecule.stereo_groups() == {}
    assert molecule.data_sgroups('COMMENT'), 'the record itself is still stored'


def test_atoms_is_a_tuple_because_one_label_covers_one_atom_or_two():
    molecule = _read(((2, 4), '*RS'), lines=_DIOL)[0]
    assert stereo_labels(molecule) == [((2, 4), (), 'RS', True)]


def test_the_container_method_is_this_function():
    """Registered onto the sealed core container by `chython.formats`, so a caller holding a molecule
    does not import the format package to read its own labels."""
    molecule = _read(((2,), '(R)'))[0]
    assert molecule.stereo_labels() == stereo_labels(molecule)


def test_a_molecule_with_no_data_group_has_no_labels():
    molecule = _read()[0]
    assert stereo_labels(molecule) == []
    assert promote_stereo_labels(molecule) is False


# ------------------------------------------------------- promote_stereo_labels(): what each word does

@mark.parametrize('datum,kind,group', [
    ('(R)', STEREO_ABS, 0),            # absolute: the record names one configuration
    ('(S)', STEREO_ABS, 0),
    ('*R', STEREO_AND, 1),             # the relative marker is what makes it a racemate
    ('*S', STEREO_AND, 1),
    ('RS', STEREO_AND, 1),
    ('*RS', STEREO_AND, 1),
])
def test_a_word_that_names_a_collection_lands_as_one(datum, kind, group):
    molecule = _read(((2,), datum))[0]
    assert molecule.stereo_groups() == {(kind, group): [2]}
    assert _rules(molecule) == ['sgroup:stereo-label-promoted']
    assert molecule.log[0].severity == REPAIRED
    assert molecule.log[0].atoms == (2,)


@mark.parametrize('datum', ['E', 'Z', '*E', '&', 'RR', 'UNDEF', ''])
def test_a_word_that_names_no_collection_is_kept_as_the_data_it_is(datum):
    """``E``/``Z`` describe an axis and no record states one beside a collection to settle the reading
    against; ``&`` occurs and is not a configuration at all.  The label survives as a data record."""
    molecule = _read(((2,), datum))[0]
    assert molecule.stereo_groups() == {}
    assert _rules(molecule) == ['sgroup:stereo-label-not-promoted']
    assert molecule.log[0].severity is LOST
    assert molecule.data_sgroups(STEREOLABEL), 'nothing was dropped -- it was not promoted'


def test_one_group_per_s_group_so_two_atoms_in_one_label_are_one_racemate():
    molecule = _read(((2, 4), '*RS'), lines=_DIOL)[0]
    assert molecule.stereo_groups() == {(STEREO_AND, 1): [2, 4]}


def test_two_labels_are_two_collections_and_the_ids_do_not_collide():
    """Two records each saying "unknown at this centre" say nothing about the other centre, so they
    cannot share a group id -- one AND group holding both would claim they are correlated."""
    molecule = _read(((2,), '*R'), ((4,), '*S'), lines=_DIOL)[0]
    assert molecule.stereo_groups() == {(STEREO_AND, 1): [2], (STEREO_AND, 2): [4]}
    assert _rules(molecule) == ['sgroup:stereo-label-promoted'] * 2


def test_an_abs_label_and_a_relative_one_coexist():
    molecule = _read(((2,), '(R)'), ((4,), '*S'), lines=_DIOL)[0]
    assert molecule.stereo_groups() == {(STEREO_ABS, 0): [2], (STEREO_AND, 1): [4]}


def test_a_label_on_an_atom_no_configuration_reaches_is_not_promoted():
    """A PROMOTION NEVER INVENTS A PARITY.  Atom 3 is the hydroxyl oxygen: it anchors nothing, so a
    membership there would name a configuration that does not exist -- the state `clean_stereo()`
    exists to remove."""
    molecule = _read(((3,), '*R'), lines=_DIOL)[0]
    assert molecule.stereo_groups() == {}
    assert _rules(molecule) == ['sgroup:stereo-label-not-promoted']
    assert molecule.log[0].severity is LOST
    assert 'does not exist' in molecule.log[0].message


def test_only_the_labelled_atoms_that_carry_a_configuration_join():
    """A pair label over one configured atom and one bare one lands on the configured atom alone,
    for the same reason: the other one has nothing to be uncertain about."""
    molecule = _read(((2, 3), '*RS'), lines=_DIOL)[0]
    assert molecule.stereo_groups() == {(STEREO_AND, 1): [2]}


def test_a_flat_drawing_promotes_nothing_at_all():
    flat = [x for x in _ALANINE]
    assert flat[10] == '  2  1  1  6  0  0  0', 'the wedge line moved in the fixture'
    flat[10] = '  2  1  1  0  0  0  0'                  # the same bond, redrawn plain
    molecule = _read(((2,), '(R)'), lines=flat)[0]
    assert molecule.parity_of(2) == 0, 'the fixture must really be flat, or this tests nothing'
    assert molecule.stereo_groups() == {}
    assert _rules(molecule) == ['sgroup:stereo-label-not-promoted']


def test_ignore_stereo_leaves_the_label_as_data():
    """Under `ignore_stereo` there is no parity for a label to belong to, so there is nothing to
    promote and the reader does not try."""
    molecule = _read(((2,), '(R)'), ignore_stereo=True)[0]
    assert molecule.parity_of(2) == 0
    assert molecule.stereo_groups() == {}
    assert _rules(molecule) == []
    assert stereo_labels(molecule) == [((2,), (), 'R', False)]


# ------------------------------------------------------------------------- a collection wins

def test_a_collection_the_record_states_is_the_partition():
    """The one rule that decides which records promotion reaches: where the file spells a collection
    the file has already said what the partition is, and where the two disagree the spec-defined
    statement is the one to believe."""
    lines = _ALANINE_V3000 + ['M  V30 BEGIN COLLECTION',
                              'M  V30 MDLV30/STERAC1 ATOMS=(1 2)',
                              'M  V30 END COLLECTION',
                              'M  V30 BEGIN SGROUP',
                              'M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME=STEREOLABEL FIELDDATA=(R)',
                              'M  V30 END SGROUP',
                              'M  V30 END CTAB', 'M  END']
    molecule = parse_v3000(lines, []).build()[0]
    assert molecule.stereo_groups() == {(STEREO_AND, 1): [2]}, 'STERAC1, not the label ABS'
    assert _rules(molecule) == ['sgroup:stereo-label-not-promoted']
    assert molecule.log[0].severity is INFO
    assert 'which is the partition' in molecule.log[0].message


def test_the_same_record_without_the_collection_block_promotes():
    """The other half of the pair, so the test above is about the collection and not about V3000."""
    lines = _ALANINE_V3000 + ['M  V30 BEGIN SGROUP',
                              'M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME=STEREOLABEL FIELDDATA=(R)',
                              'M  V30 END SGROUP',
                              'M  V30 END CTAB', 'M  END']
    molecule = parse_v3000(lines, []).build()[0]
    assert molecule.stereo_groups() == {(STEREO_ABS, 0): [2]}
    assert _rules(molecule) == ['sgroup:stereo-label-promoted']


def test_a_promoted_id_does_not_collide_with_one_the_record_already_used():
    """`STERAC1` is taken on the other centre, so the promotion takes the next free number rather than
    merging two statements the record made separately."""
    lines = _DIOL + _dat(1, (2,), '*RS') + ['M  END']
    molecule = parse_v2000(lines, []).build()[0]
    assert molecule.stereo_groups() == {(STEREO_AND, 1): [2]}
    # and again with the file's own group in the way
    v3000 = _ALANINE_V3000 + ['M  V30 BEGIN COLLECTION',
                              'M  V30 MDLV30/STERAC1 ATOMS=(1 4)',
                              'M  V30 END COLLECTION',
                              'M  V30 BEGIN SGROUP',
                              'M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME=STEREOLABEL FIELDDATA=*RS',
                              'M  V30 END SGROUP',
                              'M  V30 END CTAB', 'M  END']
    molecule = parse_v3000(v3000, []).build()[0]
    assert molecule.stereo_groups() == {(STEREO_AND, 1): [4], (STEREO_AND, 2): [2]}


# --------------------------------------------------------------------------- the record survives

def test_the_label_is_still_a_data_record_after_it_was_promoted():
    """Promotion reads the record; it does not consume it.  A V2000 round trip therefore carries the
    label back out -- and since V2000 cannot spell a collection, the label is the only carrier the
    format has and the second reading promotes it again to the same group."""
    molecule, store, _ = _read(((2,), '*RS'))
    assert molecule.stereo_groups() == {(STEREO_AND, 1): [2]}
    lines, log = emit_v2000(molecule, store)
    again = parse_v2000(lines, []).build()[0]
    assert stereo_labels(again) == [((2,), (), 'RS', True)]
    assert again.stereo_groups() == {(STEREO_AND, 1): [2]}
    assert _rules(again) == ['sgroup:stereo-label-promoted']


def test_promotion_is_idempotent_because_its_own_output_is_a_stated_collection():
    molecule = _read(((2,), '*RS'))[0]
    before = molecule.stereo_groups()
    assert promote_stereo_labels(molecule) is False
    assert molecule.stereo_groups() == before
