# -*- coding: utf-8 -*-
import pytest

try:
    from chython import inchi
    HAS_INCHI = True
except ImportError:
    HAS_INCHI = False


@pytest.mark.skipif(not HAS_INCHI, reason='libinchi not available')
def test_inchi_parse_ethanol():
    mol = inchi('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
    assert len(mol) == 3
    assert 'O' in {a.atomic_symbol for _, a in mol.atoms()}


@pytest.mark.skipif(not HAS_INCHI, reason='libinchi not available')
def test_inchi_parse_benzene():
    mol = inchi('InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H')
    assert len(mol) == 6


@pytest.mark.skipif(not HAS_INCHI, reason='libinchi not available')
def test_inchi_parse_charged():
    mol = inchi('InChI=1S/Na.H2O/h;1H2/q+1;/p-1')
    assert mol is not None
