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
from pytest import approx, mark
from chython.chemistry import crippen_logp, crippen_mr, explicify_hydrogens
from chython.chemistry._crippen import crippen_contributions
from chython.core import read_smiles


LOGP = [
    # methane: C1 0.1441 + 4 x H1 0.1230
    ('C', 0.6361),
    # ethane: 2 x C1 0.1441 + 6 x H1 0.1230
    ('CC', 1.0262),
    # benzene: 6 x C18 0.1581 + 6 x H1 0.1230
    ('c1ccccc1', 1.6866),
    # toluene: 5 x C18 0.1581 + C21 0.1360 + C8 0.08452 + 8 x H1 0.1230
    #   = 0.7905 + 0.1360 + 0.08452 + 0.9840 = 1.99502
    # Reachable only because C8 precedes C1 in crippen.tsv: in published-numbering order C1
    # `[C;z1;x0;!z4]` claims the methyl first and the sum is 2.0546 instead.  A failure here means the
    # table was re-sorted by type name, not that the arithmetic is wrong.
    ('Cc1ccccc1', 1.99502),
]

MR = [
    # methane: C1 2.503 + 4 x H1 1.057
    ('C', 6.731),
    # benzene: 6 x C18 3.350 + 6 x H1 1.057
    ('c1ccccc1', 26.442),
]


@mark.parametrize('smi,expected', LOGP)
def test_logp_is_the_sum_of_the_typed_contributions(smi, expected):
    assert crippen_logp(read_smiles(smi)) == approx(expected, abs=1e-4)


@mark.parametrize('smi,expected', MR)
def test_mr_is_the_sum_of_the_typed_contributions(smi, expected):
    assert crippen_mr(read_smiles(smi)) == approx(expected, abs=1e-4)


@mark.parametrize('smi', ['C', 'CC', 'CCO', 'c1ccccc1', 'CC(=O)Nc1ccc(O)cc1'])
def test_explicit_hydrogens_do_not_change_the_answer(smi):
    implicit = read_smiles(smi)
    explicit = read_smiles(smi)
    explicify_hydrogens(explicit)
    assert crippen_logp(explicit) == approx(crippen_logp(implicit), abs=1e-9)
    assert crippen_mr(explicit) == approx(crippen_mr(implicit), abs=1e-9)


def test_every_heavy_atom_gets_exactly_one_entry():
    m = read_smiles('CC(=O)Nc1ccc(O)cc1')
    parts = crippen_contributions(m)
    assert set(parts) == set(m.atom_numbers)
    assert sum(p[1] for p in parts.values()) == approx(crippen_logp(m))
    assert sum(p[2] for p in parts.values()) == approx(crippen_mr(m))


def test_no_carbon_is_ever_untyped_because_the_catch_all_catches_all_of_them():
    """`CS` is `[C:1]`, so no carbon can fall through, whatever it is bonded to.

    Asserts the type is not `'-'`, not merely that an entry exists: `crippen_contributions` writes
    `('-', 0.0, 0.0)` for an untypable atom, so a key check holds with or without the `CS` row.
    Every carbon here is one an early carbon row cannot claim -- `C[Se]C` and `C[Zn]C`'s methyls fail
    `C1`'s `x0`, and the alkyne, alkene and epoxide carbons are outside the sp3-C-and-H-only block.
    """
    for smi in ('C', 'C#C', 'C=C', 'CC(=O)O', 'C[Se]C', 'C[Zn]C', 'c1ccccc1', 'C1CO1'):
        m = read_smiles(smi)
        parts = crippen_contributions(m)
        for atom in m.atoms():
            if atom.element == 6:
                assert parts[atom.n][0] != '-', (smi, atom.n)


def test_an_element_outside_the_papers_set_is_answered_about_and_logged():
    """Gadolinium has no Wildman-Crippen type and no catch-all can reach it.

    Every catch-all is element-scoped -- `CS` is `[C:1]`, `NS` a nitrogen, `OS` an oxygen -- so an
    element outside the paper's set is typed `'-'`, contributes zero to both quantities, and puts one
    record on `m.log`: answered about, never refused.

    The element must be a lanthanide, which the paper genuinely omits.  Selenium will not do: the
    publication's last `Me1` alternative is `[#34,#52,#84]`, so Se, Te and Po are Me1 at -0.3808.
    """
    m = read_smiles('[Gd]')
    parts = crippen_contributions(m)
    log = m.log
    assert set(parts) == set(m.atom_numbers)
    gd = next(a.n for a in m.atoms() if a.element == 64)
    assert parts[gd] == ('-', 0.0, 0.0)
    assert [r.rule for r in log] == ['crippen:untyped']
    assert log[0].atoms == (gd,)


def test_selenium_is_me1_because_the_publication_says_so():
    """`Se` and `Te` are `Me1` per the publication, not untyped elements.

    Its own test because the failure it prevents is invisible: a missing element in a 93-entry list.
    """
    for smi, n in (('[Se]', 34), ('[Te]', 52)):
        parts = crippen_contributions(read_smiles(smi))
        assert parts[next(a.n for a in read_smiles(smi).atoms() if a.element == n)][0] == 'Me1'


def test_a_dative_oxide_oxygen_is_typed_rather_than_left_untyped():
    """O5 and O6 claim the anionic oxide oxygen, which is where the publication puts it.

    The measurement is O7's own pattern: "other anionic oxygen" excludes a `#7` and a `#16` neighbour,
    so the paper routes an anionic oxygen on nitrogen to O5 and one on sulfur to O6 rather than to O7 --
    and O5's transcribed probe, `C[N+](=O)[O-]`, is a nitro compound for the same reason.  With the
    mapped oxygen charge-unstated neither row reached the atom both were written for, and it fell past
    the OS catch-all, which is charge-unstated too, to `'-'`.

    A phosphorus neighbour is absent here on purpose: O7 excludes only `#7` and `#16`, so the oxygen of
    `C[P+](C)(C)[O-]` is O7 by the paper's own construction and is not this defect.
    """
    expected = {'C[N+](C)(C)[O-]': ('O5',),                       # trimethylamine N-oxide
                'c1cc[n+]([O-])cc1': ('O5',),                     # pyridine N-oxide
                'C[N+](=NC)[O-]': ('O5',),                        # azoxymethane
                'c1ccccc1[N+](=O)[O-]': ('O5', 'O5'),             # nitrobenzene: both oxygens
                'C[S+]([O-])C': ('O6',),                          # the charge-separated sulfoxide
                'c1ccccc1S(=O)(=O)[O-]': ('O6', 'O6', 'O6')}      # benzenesulfonate: all three
    for smi, types in expected.items():
        m = read_smiles(smi)
        parts = crippen_contributions(m)
        assert tuple(parts[a.n][0] for a in m.atoms() if a.element == 8) == types, smi
        assert 'crippen:untyped' not in [r.rule for r in m.log], smi


def test_the_phosphine_oxide_oxygen_stays_o7():
    """The neighbour O7 does not exclude, pinned so widening O5 and O6 cannot quietly reach it."""
    m = read_smiles('C[P+](C)(C)[O-]')
    assert crippen_contributions(m)[next(a.n for a in m.atoms() if a.element == 8)][0] == 'O7'


def test_the_container_properties_agree_with_the_functions():
    m = read_smiles('Cc1ccccc1')
    assert m.crippen_logp == approx(crippen_logp(m))
    assert m.crippen_mr == approx(crippen_mr(m))


def test_it_is_renumbering_invariant():
    assert crippen_logp(read_smiles('Cc1ccccc1')) == approx(crippen_logp(read_smiles('c1ccccc1C')))
