# -*- coding: utf-8 -*-
import pytest

from chython.core import MoleculeContainer
from chython.core import STEREO_ABS, STEREO_AND, STEREO_OR, STEREO_UNSPECIFIED
from chython.core import WEDGE_DOWN, WEDGE_EITHER, WEDGE_NONE, WEDGE_UP


def test_molecule_without_coordinates_reports_none():
    m = MoleculeContainer()
    a = m.add_atom(6)
    assert m.has_coordinates is False
    assert m.xy_of(a) is None


def test_coordinates_survive_the_apply_at_1e4_resolution():
    m = MoleculeContainer()
    a1, a2 = m.add_atom(6), m.add_atom(8)
    m.add_bond(a1, a2, 1)
    m.set_xy(a1, 0.0, 0.0)
    m.set_xy(a2, 1.2345, -0.8765)
    assert m.has_coordinates is True
    assert m.xy_of(a1) == (0.0, 0.0)
    x, y = m.xy_of(a2)
    assert x == pytest.approx(1.2345, abs=5e-5)
    assert y == pytest.approx(-0.8765, abs=5e-5)


def test_partial_coordinates_default_to_origin():
    m = MoleculeContainer()
    a1, a2 = m.add_atom(6), m.add_atom(6)
    m.add_bond(a1, a2, 1)
    m.set_xy(a1, 3.0, 4.0)
    assert m.xy_of(a1) == (3.0, 4.0)
    assert m.xy_of(a2) == (0.0, 0.0)


def test_coordinates_out_of_fixed_point_range_raise():
    m = MoleculeContainer()
    a = m.add_atom(6)
    with pytest.raises(ValueError):
        m.set_xy(a, 300000.0, 0.0)


def test_coordinates_survive_mutation():
    m = MoleculeContainer()
    a1, a2 = m.add_atom(6), m.add_atom(8)
    m.add_bond(a1, a2, 1)
    m.set_xy(a1, 1.5, 2.5)
    m.set_xy(a2, 3.5, 4.5)
    with m.edit():
        m.set_charge(a2, -1)
    assert m.xy_of(a1) == (1.5, 2.5)
    assert m.xy_of(a2) == (3.5, 4.5)


def test_coordinates_survive_compaction():
    m = MoleculeContainer()
    a = m.add_atom(6)
    b = m.add_atom(6)
    c = m.add_atom(6)
    m.set_xy(a, 1.0, 0.0)
    m.set_xy(b, 2.0, 0.0)
    m.set_xy(c, 3.0, 0.0)
    m.delete_atom(b)
    assert m.xy_of(a) == (1.0, 0.0)
    assert m.xy_of(c) == (3.0, 0.0)


def _chiral_center():
    m = MoleculeContainer()
    c = m.add_atom(6)
    f = m.add_atom(9)
    cl = m.add_atom(17)
    br = m.add_atom(35)
    for x in (f, cl, br):
        m.add_bond(c, x, 1)
    return m, c, f, cl, br


def test_wedge_is_directed_and_the_twin_stays_clean():
    m, c, f, cl, br = _chiral_center()
    m.set_wedge(c, f, WEDGE_UP)
    assert m.wedge_of(c, f) == WEDGE_UP
    assert m.wedge_of(f, c) == WEDGE_NONE


def test_wedge_configures_no_parity_on_the_narrow_atom():
    """A wedge is geometry only; it must not write a parity (Ruling F54).

    After set_wedge, parity_of must be 0 (unset) and stereo_of must be False on the narrow atom.
    The non-narrow atom (f) also stays False.
    """
    m, c, f, cl, br = _chiral_center()
    m.set_wedge(c, f, WEDGE_DOWN)
    assert m.parity_of(c) == 0
    assert m.stereo_of(c) is False
    assert m.stereo_of(f) is False


def test_multiple_wedges_are_listed_narrow_first():
    m, c, f, cl, br = _chiral_center()
    m.set_wedge(c, f, WEDGE_UP)
    m.set_wedge(c, cl, WEDGE_DOWN)
    m.set_wedge(c, br, WEDGE_EITHER)
    assert sorted(m.wedges()) == sorted([(c, f, WEDGE_UP), (c, cl, WEDGE_DOWN),
                                         (c, br, WEDGE_EITHER)])


def test_wedges_survive_mutation():
    m, c, f, cl, br = _chiral_center()
    m.set_wedge(c, f, WEDGE_UP)
    with m.edit():
        m.set_charge(f, -1)
    assert m.wedge_of(c, f) == WEDGE_UP


def test_wedge_on_a_nonexistent_bond_raises():
    m, c, f, cl, br = _chiral_center()
    with pytest.raises(KeyError):
        m.set_wedge(f, cl, WEDGE_UP)


def test_invalid_wedge_value_raises():
    m, c, f, cl, br = _chiral_center()
    with pytest.raises(ValueError):
        m.set_wedge(c, f, 4)


def test_molecule_without_wedges_reports_none_everywhere():
    m, c, f, cl, br = _chiral_center()
    assert m.wedges() == []
    assert m.wedge_of(c, f) == WEDGE_NONE


def test_wedge_is_dropped_when_its_bond_is_deleted():
    m = MoleculeContainer()
    a = m.add_atom(6)
    b = m.add_atom(6)
    m.add_bond(a, b, 1)
    m.set_wedge(a, b, WEDGE_UP)
    m.delete_bond(a, b)
    assert m.wedges() == []
    assert a in m.atom_numbers
    assert b in m.atom_numbers


def test_setting_a_wedge_then_deleting_its_bond_in_one_scope_raises():
    # A contradictory edit: wedge and delete_bond in one scope. The raise is deliberate.
    m = MoleculeContainer()
    a = m.add_atom(6)
    b = m.add_atom(6)
    m.add_bond(a, b, 1)
    with pytest.raises(KeyError):
        with m.edit():
            m.set_wedge(a, b, WEDGE_UP)
            m.delete_bond(a, b)


def _two_centres():
    m = MoleculeContainer()
    ids = [m.add_atom(6) for _ in range(4)]
    m.add_bond(ids[0], ids[1], 1)
    m.add_bond(ids[1], ids[2], 1)
    m.add_bond(ids[2], ids[3], 1)
    return m, ids


def test_no_stereo_groups_means_unspecified_everywhere():
    m, ids = _two_centres()
    assert m.has_stereo_groups is False
    assert m.stereo_group_of(ids[0]) == (STEREO_UNSPECIFIED, 0)
    assert m.stereo_groups() == {}


def test_abs_ignores_the_group_id():
    m, ids = _two_centres()
    m.set_stereo_group(ids[1], STEREO_ABS)
    assert m.has_stereo_groups is True
    assert m.stereo_group_of(ids[1]) == (STEREO_ABS, 0)


def test_or_and_and_keep_distinct_group_ids():
    m, ids = _two_centres()
    m.set_stereo_group(ids[1], STEREO_OR, 1)
    m.set_stereo_group(ids[2], STEREO_AND, 1)
    assert m.stereo_group_of(ids[1]) == (STEREO_OR, 1)
    assert m.stereo_group_of(ids[2]) == (STEREO_AND, 1)
    # same group number, different kind => different collections
    assert m.stereo_groups() == {(STEREO_OR, 1): [ids[1]], (STEREO_AND, 1): [ids[2]]}


def test_members_of_one_group_are_collected_together():
    m, ids = _two_centres()
    m.set_stereo_group(ids[1], STEREO_AND, 2)
    m.set_stereo_group(ids[2], STEREO_AND, 2)
    assert m.stereo_groups() == {(STEREO_AND, 2): [ids[1], ids[2]]}


def test_group_id_bounds():
    m, ids = _two_centres()
    m.set_stereo_group(ids[1], STEREO_OR, 63)
    with pytest.raises(ValueError):
        m.set_stereo_group(ids[2], STEREO_OR, 64)
    with pytest.raises(ValueError):
        m.set_stereo_group(ids[2], STEREO_OR, 0)     # OR/AND need a real group
    with pytest.raises(ValueError):
        m.set_stereo_group(ids[2], 4, 1)             # unknown kind


def test_stereo_groups_survive_mutation():
    m, ids = _two_centres()
    m.set_stereo_group(ids[1], STEREO_AND, 3)
    with m.edit():
        m.set_charge(ids[0], 1)
    assert m.stereo_group_of(ids[1]) == (STEREO_AND, 3)
