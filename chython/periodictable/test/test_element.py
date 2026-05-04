# -*- coding: utf-8 -*-
import pytest
from chython.periodictable import C, N, O, H, S, F, Cl, Br
from chython.periodictable.base.element import Element, _compiled_valence_rules, _compiled_saturation_rules


def test_element_from_symbol():
    cls = Element.from_symbol('C')
    assert cls is C


def test_element_from_atomic_number():
    cls = Element.from_atomic_number(6)
    assert cls is C
    cls = Element.from_atomic_number(7)
    assert cls is N


def test_element_from_atomic_number_invalid():
    with pytest.raises(ValueError):
        Element.from_atomic_number(999)


def test_carbon_valence_rules():
    c = C()
    rules = c.valence_rules(4)
    assert len(rules) > 0
    # Carbon with 4 bonds should have 0 implicit H
    assert any(h == 0 for _, _, h in rules)


def test_carbon_valence_rules_partial():
    c = C()
    rules = c.valence_rules(3)
    # Carbon with 3 bonds should have 1 implicit H
    assert any(h == 1 for _, _, h in rules)


def test_nitrogen_charged_valence():
    n = N(charge=1)
    rules = n.valence_rules(4)
    assert len(rules) > 0


def test_invalid_valence():
    from chython.exceptions import ValenceError
    c = C()
    with pytest.raises(ValenceError):
        c.valence_rules(99)


def test_saturation_rules():
    rules = _compiled_saturation_rules(C)
    assert len(rules) > 0
    # Common valence 4 should be there
    assert any(v == 4 for _, _, v, _, _ in rules)


def test_element_properties():
    c = C()
    assert c.atomic_symbol == 'C'
    assert c.atomic_number == 6
    assert c.charge == 0
    assert c.is_radical is False


def test_element_charge():
    n = N(charge=-1)
    assert n.charge == -1


def test_element_isotope():
    c = C(isotope=13)
    assert c.isotope == 13


def test_element_invalid_isotope():
    with pytest.raises(ValueError):
        C(isotope=999)


def test_element_copy():
    c = C(isotope=13, charge=1)
    c2 = c.copy()
    assert c2.isotope == 13
    assert c2.charge == 1
    assert c2 is not c


def test_compiled_valence_rules_cached():
    # Call twice, should return same object (cached)
    r1 = _compiled_valence_rules(C)
    r2 = _compiled_valence_rules(C)
    assert r1 is r2


def test_different_elements_different_rules():
    r_c = _compiled_valence_rules(C)
    r_n = _compiled_valence_rules(N)
    assert r_c is not r_n
