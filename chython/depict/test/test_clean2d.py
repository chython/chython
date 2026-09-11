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
"""Gates on the default 2D layout backend, `clean2d_engine = 'smilesdrawer'`.

A layout is unchanged by rotation, reflection and translation and the bundle picks arbitrarily between
symmetric alternatives, so literal coordinates are not a property of the algorithm.  What is asserted is
that every molecule gets a layout, that no two atoms share a point, and that rings come out even.
"""
from math import dist, hypot
from re import findall
from pytest import approx, mark, raises, skip
from chython.core import ReactionContainer, read_smiles as smiles
from chython.depict.layout import molecule as layout_molecule
from chython.depict._config import get_clean2d_engine
from chython.depict.layout.molecule import ctx, _clean2d_tree
from chython.exceptions import ImplementationError


# `chython.core` and not `chython`: nothing under `depict/` reads a name off the facade.  Importing
# `chython.depict` for `ctx` is what registers `clean2d` onto the container.


def test_the_default_engine_is_loaded():
    """`quickjs-ng` is a required dependency, so a context that failed to build is a failure here.

    This is a diagnostic and not a safety net: with `quickjs` made unimportable, 123 of this file's
    tests already fail and only 10 skip, so the state was never going to pass unnoticed.  What it buys
    is the first failure in file order naming the cause -- `layout/molecule.py` catches a bare
    `Exception` around the context build so that an install without the engine still imports, which
    means the other 123 report a missing layout rather than a missing engine.
    """
    assert get_clean2d_engine() == 'smilesdrawer', 'the default engine is the shipped one'
    assert ctx is not None, ('quickjs did not load, so `clean2d_engine = \'smilesdrawer\'` raises: '
                             'either quickjs-ng is not installed or it is built against an '
                             'incompatible libquickjs')
    # A loaded context is not a working one -- it is built from `clean2d.js` plus a shim, and a bundle
    # that parsed can still fail on the first call.
    mol = smiles('c1ccccc1')
    mol.clean2d()
    assert len(mol.coordinates()) == len(mol), 'the engine answered a point for every atom'


# (name, smiles) -- public structures only.
flat = [
    ('benzene', 'c1ccccc1'),
    ('naphthalene', 'c1ccc2ccccc2c1'),
    ('anthracene', 'c1ccc2cc3ccccc3cc2c1'),
    ('biphenyl', 'c1ccc(cc1)-c1ccccc1'),
    ('pyridine', 'c1ccncc1'),
    ('indole', 'c1ccc2[nH]ccc2c1'),
    ('quinoxaline', 'c1ccc2nccnc2c1'),
    ('carbazole', 'c1ccc2[nH]c3ccccc3c2c1'),
    ('cyclohexane', 'C1CCCCC1'),
    ('cyclopropane', 'C1CC1'),
    ('spiro[4.5]decane', 'C1CCC2(CC1)CCCC2'),
    ('spirooxindole', 'O=C1Nc2ccccc2C11CCNCC1'),
    ('cyclohexanone_ethylene_ketal', 'O=C1CCC2(CC1)OCCO2'),
    ('aspirin', 'CC(=O)Oc1ccccc1C(=O)O'),
    ('caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'),
    ('ibuprofen', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'),
    ('naproxen', 'COc1ccc2cc(ccc2c1)C(C)C(=O)O'),
    ('nicotine', 'CN1CCCC1c1cccnc1'),
    ('atenolol', 'CC(C)NCC(O)COc1ccc(CC(N)=O)cc1'),
    ('penicillin_g', 'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O'),
    ('ampicillin', 'CC1(C)SC2C(NC(=O)C(N)c3ccccc3)C(=O)N2C1C(=O)O'),
    ('chlorpromazine', 'CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21'),
    ('glucose', 'OCC1OC(O)C(O)C(O)C1O'),
    ('maltose', 'OCC1OC(OC2C(O)C(O)C(O)OC2CO)C(O)C(O)C1O'),
    ('cholesterol', 'CC(C)CCCC(C)C1CCC2(C)C1CCC1C2CC=C2CC(O)CCC12C'),
    ('progesterone_core', 'CC(=O)C1CCC2(C)C1CCC1C2CCC2(C)C1CCC2=O'),
    ('macrolactone', 'CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)CC(C)C(=O)C(C)C(O)C1C'),
    ('porphine', 'c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3'),
    ('tryptophan', 'NC(Cc1c[nH]c2ccccc12)C(=O)O'),
    ('quinuclidine', 'C1CN2CCC1CC2'),
    ('anthraquinone_dimethoxy', 'COc1cc2c(cc1OC)C(=O)c1ccccc1C2=O'),
    ('stilbene', 'c1cc(ccc1)C=Cc1ccccc1'),
    ('boc_piperazine', 'CC(C)(C)OC(=O)N1CCNCC1'),
    ('sulfonamide', 'FC(F)(F)c1ccc(cc1)S(=O)(=O)N'),
]

# Cage systems.  A polycyclic cage has no faithful planar embedding, so `smilesdrawer` legitimately
# superimposes atoms on these.  They are still exercised -- the engine must return a layout without
# raising -- but they are exempt from the separation and evenness gates.
cages = [
    ('adamantane', 'C1C2CC3CC1CC(C2)C3'),
    ('cubane', 'C12C3C4C5C(C14)C2CC35'),
    ('norbornane', 'C1CC2CCC1CC2'),
]


def _bond_lengths(mol, plane=None):
    """Every bond's length, read off the molecule's STORED coordinates, or off a given plane."""
    if plane is None:
        plane = mol.coordinates()
    lengths = []
    for bond in mol.bonds():
        nx, ny = plane[bond.n]
        mx, my = plane[bond.m]
        lengths.append(hypot(nx - mx, ny - my))
    return lengths


def _min_separation(mol):
    xy = list(mol.coordinates().values())
    return min(dist(xy[i], xy[j]) for i in range(len(xy)) for j in range(i + 1, len(xy)))


@mark.parametrize('name,smi', flat + cages, ids=[n for n, _ in flat + cages])
def test_layout_is_produced(name, smi):
    """every molecule gets finite, non-degenerate coordinates"""
    mol = smiles(smi)
    mol.clean2d(engine='smilesdrawer')

    xy = list(mol.coordinates().values())
    assert len(xy) == len(mol)
    for x, y in xy:
        assert x == x and y == y, f'{name}: NaN coordinate'
        assert abs(x) < 1e6 and abs(y) < 1e6, f'{name}: coordinate ran away'
    assert max(hypot(x, y) for x, y in xy) > 0., f'{name}: whole molecule collapsed to the origin'


@mark.parametrize('name,smi', flat, ids=[n for n, _ in flat])
def test_no_overlapping_atoms(name, smi):
    """no two atoms closer than half a bond -- the layout is readable"""
    mol = smiles(smi)
    mol.clean2d(engine='smilesdrawer')

    lengths = _bond_lengths(mol)
    mean = sum(lengths) / len(lengths)
    assert _min_separation(mol) / mean > .45, f'{name}: atoms collide'


@mark.parametrize('name,smi', flat, ids=[n for n, _ in flat])
def test_bonds_are_even(name, smi):
    """bond lengths stay close to uniform -- rings are not distorted"""
    mol = smiles(smi)
    mol.clean2d(engine='smilesdrawer')

    lengths = _bond_lengths(mol)
    mean = sum(lengths) / len(lengths)
    spread = (sum((v - mean) ** 2 for v in lengths) / len(lengths)) ** .5 / mean
    assert spread < .2, f'{name}: bond lengths vary by {spread:.0%} of the mean'


def test_rescaled_to_standard_bond_length():
    """`clean2d` normalises to chython's 0.825 mean bond length whatever the engine returns

    `abs=1e-9` on the RETURNED plane, `2e-4` on the stored one: the arena keeps a coordinate as `xy_t`,
    an `int32_t` scaled by 10000, so storing quantises onto a 1e-4 grid, which bounds the mean bond
    length's drift at 1.5e-4 (measured 1.1e-5).  Tightening the second would assert `xy_t` is a double.
    """
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    plane = mol.layout2d(engine='smilesdrawer')

    lengths = _bond_lengths(mol, plane)
    assert sum(lengths) / len(lengths) == approx(.825, abs=1e-9)

    mol.clean2d(engine='smilesdrawer')
    stored = _bond_lengths(mol)
    assert sum(stored) / len(stored) == approx(.825, abs=2e-4)


# `rescale2d` is that normalisation on its own, applied to coordinates the molecule already carries --
# the operation a plane from a drawing editor needs, where a redraw would throw the drawing away.


def test_rescale2d_normalises_a_plane_that_came_in_at_another_scale():
    """The whole point: an editor's bond length becomes chython's 0.825 and no atom moves relative to
    any other.  `abs=2e-4` for the arena's 1e-4 coordinate grid, as above."""
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    mol.clean2d(engine='smilesdrawer')
    with mol.edit():
        for n, (x, y) in mol.coordinates().items():
            mol.set_xy(n, x * 4., y * 4.)

    assert mol.rescale2d()
    stored = _bond_lengths(mol)
    assert sum(stored) / len(stored) == approx(.825, abs=2e-4)


def test_rescale2d_answers_false_and_stores_nothing_when_there_is_no_plane_to_rescale():
    """Two ways to have no scale, and both are the same answer.

    A molecule with no coordinates reads as every atom at the origin, so its bonds have no length; a
    single atom has no bonds at all.  Neither is rescalable and neither is an error -- and the first
    must not leave a stored plane of zeros behind, which is what would make `has_layout` lie.
    """
    mol = smiles('CCO')
    assert not mol.rescale2d()
    assert not mol.has_layout

    lone = smiles('[Na+]')
    assert not lone.rescale2d()
    assert not lone.has_layout


def test_disconnected_components_are_separated():
    """`.`-joined components are laid out side by side, not on top of each other"""
    mol = smiles('CCO.c1ccccc1')
    mol.clean2d(engine='smilesdrawer')

    assert mol.connected_components_count == 2
    a, b = mol.connected_components
    xy = mol.coordinates()
    assert max(xy[n][0] for n in a) < min(xy[n][0] for n in b)


def test_layout_is_deterministic():
    """the same molecule laid out twice gives the same coordinates

    The engine is one long-lived QuickJS context, so this also gates that a layout leaves no state behind.
    """
    smi = 'CC(C)CCCC(C)C1CCC2(C)C1CCC1C2CC=C2CC(O)CCC12C'
    first = smiles(smi)
    first.clean2d(engine='smilesdrawer')
    expected = first.coordinates()
    for _ in range(3):
        again = smiles(smi)
        again.clean2d(engine='smilesdrawer')
        got = again.coordinates()
        assert got.keys() == expected.keys()
        for n, xy in expected.items():
            assert got[n] == approx(xy, abs=1e-9)


def test_engine_survives_many_layouts():
    """the shared context does not accumulate heap across calls

    `clean2d` skips the binding's per-call collection and relies on QuickJS's own threshold collector;
    were that to stop holding, the default engine would leak once per drawn structure.
    """
    if ctx is None:
        skip('quickjs is not installed')

    mol = smiles('c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3')
    tree, _ = _clean2d_tree(mol)

    ctx(tree, run_gc=False)
    ctx.gc()
    baseline = ctx.memory()['memory_used_size']
    for _ in range(200):
        ctx(tree, run_gc=False)
    assert ctx.memory()['memory_used_size'] < baseline + 4 * 1024 * 1024


def test_unknown_engine_rejected():
    mol = smiles('c1ccccc1')
    with raises(ValueError):
        mol.clean2d(engine='no-such-engine')


def test_depict_must_not_mutate_coordinates():
    """a renderer computes its layout into a temporary: asking for a picture is not asking for geometry

    A caller that wants coordinates calls `clean2d()`.
    """
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    before = mol.coordinates()
    assert not before, 'expected a molecule with no coordinates'

    svg = mol.depict()
    assert '<path' in svg, 'nothing was drawn, so the test proves nothing'

    assert mol.coordinates() == before, 'depict() rewrote the coordinates it was given'


def test_repr_svg_must_not_mutate_coordinates():
    """the Jupyter path is the one that fires without anyone asking to draw"""
    mol = smiles('c1ccc2[nH]ccc2c1')
    before = mol.coordinates()

    mol._repr_svg_()

    assert mol.coordinates() == before


def test_reaction_depict_must_not_mutate_coordinates():
    """a reaction is arranged left to right for drawing; that arrangement is not the reaction

    The arrow and the signs are RETURNED and a reaction has nowhere to keep them, so what is observable
    is whether the molecules' coordinates survive.
    """
    rxn = smiles('CC(=O)O.NCC>>CC(=O)NCC')
    before = [m.coordinates() for m in rxn.molecules()]

    svg = rxn.depict()
    assert '<path' in svg

    assert [m.coordinates() for m in rxn.molecules()] == before


def test_wedges_survive_drawing_without_stored_coordinates():
    """which bond carries a wedge is a property of the drawing, not of the molecule

    Wedge direction has to be derived from the positions being RENDERED: read off unset stored
    coordinates it sees every atom at the origin, the pyramid sign comes out zero and the marks vanish.
    """
    mol = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O')     # glucopyranose
    assert not mol.coordinates()

    unlaid = _filled_paths(mol.depict())
    assert unlaid, 'no wedge was drawn for a molecule with five stereocentres'

    # and the same molecule laid out first draws the same number of them
    mol.clean2d()
    assert _filled_paths(mol.depict()) == unlaid


def _filled_paths(svg):
    """How many wedges the picture holds: a `<path>` with a solid fill, counted off the markup.

    `<path>` and not `fill="#000000"` alone -- a black atom label is a `<text>` with the same fill.
    """
    return len([tag for tag in findall(r'<path[^>]*>', svg) if 'fill="#000000"' in tag])


def test_clean2d_still_keeps_the_layout():
    """the opposite of the above: an explicit `clean2d()` does store, and the arrow is RETURNED"""
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    mol.clean2d()
    assert mol.has_layout
    assert len(mol.coordinates()) == len(mol)

    rxn = smiles('CC(=O)O.NCC>>CC(=O)NCC')
    arrow, signs = rxn.clean2d()
    assert arrow is not None
    for m in rxn.molecules():
        assert m.has_layout, 'a member of the reaction was left without coordinates'


def _assert_signs_sit_in_the_gaps(rxn, planes, signs):
    """Every `+` lies strictly between the two members it separates, on both sides of the arrow.

    Positional rather than a count, which cannot see a sign between the wrong pair, and read off the
    SHIFTED planes so the check is on where the `+` lands among the atoms.
    """
    reactants, agents, products = len(rxn.reactants), len(rxn.agents), len(rxn.products)
    spans = [(min(x for x, _ in p.values()), max(x for x, _ in p.values())) for p in planes]

    gaps = []
    for i in range(reactants - 1):
        gaps.append((spans[i][1], spans[i + 1][0]))
    for i in range(reactants + agents, reactants + agents + products - 1):
        gaps.append((spans[i][1], spans[i + 1][0]))

    assert len(signs) == len(gaps), 'one sign per gap between two members of a side, and no others'
    for (x, y), (left, right) in zip(signs, gaps):
        assert left < x < right, f'a + at {x} is not in the gap ({left}, {right}) it separates'
        assert y == 0., 'the row is the y = 0 axis, so a sign sits on it'


def test_a_reaction_layout_places_an_arrow_between_the_sides():
    rxn = ReactionContainer([smiles('CCO'), smiles('CC(=O)O')], [smiles('CCOC(C)=O')])
    planes, arrow, signs = rxn.layout2d()
    x1, x2, y = arrow
    assert x2 > x1, 'the arrow spans left to right'
    _assert_signs_sit_in_the_gaps(rxn, planes, signs)


def test_a_multi_product_reaction_signs_its_product_side_too():
    """the products loop's `if amount:` body, which nothing else in the suite enters

    Ester hydrolysis is two members on each side, so one `+` is expected on each and the reactant-side
    count cannot stand in for both.
    """
    rxn = smiles('CC(=O)OCC.O>>CC(=O)O.CCO')
    planes, arrow, signs = rxn.layout2d()
    arrow_min, arrow_max, _ = arrow

    assert len(rxn.reactants) == 2 and len(rxn.products) == 2
    assert len(signs) == 2, 'one + between the reactants and one between the products'
    _assert_signs_sit_in_the_gaps(rxn, planes, signs)

    reactant_sign, product_sign = signs
    assert reactant_sign[0] < arrow_min, 'the reactant + is left of the arrow'
    assert product_sign[0] > arrow_max, 'the product + is right of the arrow'


def test_an_agent_sits_inside_the_arrow_span():
    """the agent branch, and `_shift_plane_min`, which nothing else in the suite enters

    An agent is drawn ON the arrow, so `_position` insets it by `+.4` at each end and floors the arrow's
    advance so the arrow is never shorter than what it carries.
    """
    rxn = smiles('CC(=O)O.NCC>[Pd]>CC(=O)NCC')
    planes, arrow, signs = rxn.layout2d()
    arrow_min, arrow_max, _ = arrow

    assert len(rxn.agents) == 1
    agent_plane = planes[len(rxn.reactants)]
    for x, y in agent_plane.values():
        assert arrow_min < x < arrow_max, 'the agent is drawn outside the arrow it sits on'

    assert arrow_max - arrow_min >= 2., 'the arrow is shorter than its minimum span'
    _assert_signs_sit_in_the_gaps(rxn, planes, signs)


def test_a_reactions_second_layout_is_stable():
    """what `reaction.clean2d`'s docstring claims: one quantisation step, then exact

    The first call arranges unquantised planes and the second the stored, rounded ones, so the arrow
    moves once by under a grid step and never again.
    """
    rxn = smiles('CC(=O)O.NCC>[Pd]>CC(=O)NCC')
    first, _ = rxn.clean2d()
    stored = [m.coordinates() for m in rxn.molecules()]

    second, _ = rxn.clean2d()
    assert [m.coordinates() for m in rxn.molecules()] == stored, 'the stored planes moved'
    assert second == approx(first, abs=1e-4), 'the arrow moved by more than one quantisation step'

    assert rxn.clean2d()[0] == second, 'the arrow is still moving on the third call'


def test_a_single_atom_plane_is_floats():
    """a plane is `{n: (float, float)}` for every atom, including one laid out at the origin

    QuickJS returns an integral JS number as a Python `int`, and a one-atom molecule lands on an
    integral origin, so without the cast one atom in a plane is a pair of `int`s.
    """
    for smi in ('[Na+]', 'CCO.[Na+]', 'c1ccccc1'):
        plane = smiles(smi).layout2d(engine='smilesdrawer')
        for xy in plane.values():
            assert [type(v) for v in xy] == [float, float], f'{smi}: {xy!r} is not a pair of floats'


def _stretched(mol):
    """Store a plane no engine would produce, and return it.

    Non-degenerate on purpose: `has_layout` is a bounding-box test, so a plane a `force=False` call must
    KEEP has to read as a layout.  An integer lattice survives the arena's 1e-4 grid exactly.
    """
    plane = {}
    ids = list(mol)
    for i, n in enumerate(ids):
        plane[n] = (float(i), float(i % 3))
    with mol.edit():
        for n, (x, y) in plane.items():
            mol.set_xy(n, x, y)
    return plane


def test_clean2d_without_force_keeps_an_existing_layout():
    """`has_layout` true, `force` false -- a no-op, and the correct answer to "make sure it has one" """
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    stretched = _stretched(mol)
    assert mol.has_layout

    mol.clean2d(engine='smilesdrawer')
    assert mol.coordinates() == stretched, 'clean2d() overwrote a layout it was not asked to replace'


def test_clean2d_with_force_replaces_an_existing_layout():
    """`has_layout` true, `force` true -- replacing a layout is opt-in"""
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    mol.clean2d(engine='smilesdrawer')
    laid = mol.coordinates()
    _stretched(mol)

    mol.clean2d(engine='smilesdrawer', force=True)
    assert mol.coordinates() == laid, 'force=True did not relay the molecule'


def test_clean2d_lays_out_a_molecule_with_no_layout():
    """`has_layout` false, `force` false -- the compute branch, which is the common case"""
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    assert not mol.has_layout

    mol.clean2d(engine='smilesdrawer')
    assert mol.has_layout
    lengths = _bond_lengths(mol)
    assert sum(lengths) / len(lengths) == approx(.825, abs=2e-4)


def test_layout2d_returns_the_stored_plane_unless_forced():
    """the same question on the non-storing entry point: from storage unless `force=True`"""
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    computed = mol.layout2d(engine='smilesdrawer')
    assert not mol.has_layout, 'layout2d() stored its result'

    stretched = _stretched(mol)
    assert mol.layout2d(engine='smilesdrawer') == stretched

    relaid = mol.layout2d(engine='smilesdrawer', force=True)
    assert relaid.keys() == computed.keys()
    for n, xy in computed.items():
        assert relaid[n] == approx(xy, abs=1e-9)
    assert mol.coordinates() == stretched, 'layout2d(force=True) stored its result'


def test_a_partial_layout_registration_is_refused_at_the_hook():
    """the five layout functions are one registration, and the setter is where that is enforced

    Only the REFUSAL is exercised: a call that got past it would replace this process's real
    registration with the stub, and every layout test after it would measure that instead.
    """
    from chython.core._core import _set_depict_fns

    def stub(*args, **kwargs):
        raise AssertionError('the stub must never be reachable: the call above has to refuse first')

    layout_names = ('clean2d', 'layout2d', 'rescale2d', 'reaction_clean2d', 'reaction_layout2d')
    for offered in layout_names:
        with raises(ValueError, match='ONE registration'):
            _set_depict_fns(**{offered: stub})

    # every 4-of-5 combination: each omitted name must appear in the error message
    for omitted in layout_names:
        with raises(ValueError, match=omitted):
            _set_depict_fns(**{n: stub for n in layout_names if n != omitted})

    # and a call offering NOTHING from this group leaves it alone rather than raising -- that is what
    # lets a later group (the drawing entry points) register itself in a call of its own.
    _set_depict_fns()
    smiles('c1ccccc1').clean2d(engine='smilesdrawer')


# An EXPLICIT hydrogen is the atom the default engine will not lay out: smiles-drawer's graph builder
# gives a node an index only when the element is not H or the node is a lone root, so a written-out `[H]`
# costs a point that never comes back.  It is how a wedge-bearing structure comes out of an MDL file.


def test_an_explicit_hydrogen_gets_a_point_of_its_own():
    """the whole molecule is laid out, and NO TWO ATOMS SHARE A POINT

    The separation is the half that discriminates: a mis-assigned plane is otherwise well formed, so a
    test that only asked for coordinates passes on it.
    """
    if ctx is None:
        skip('quickjs is not installed')

    mol = smiles('F[C@]([H])(Cl)Br')
    mol.clean2d(engine='smilesdrawer')

    xy = mol.coordinates()
    assert xy.keys() == set(mol)
    lengths = _bond_lengths(mol)
    assert _min_separation(mol) / (sum(lengths) / len(lengths)) > .45, 'two atoms collided'


def test_the_deferred_hydrogen_sits_a_bond_length_from_its_neighbour():
    """placed at a PLAUSIBLE distance, asserted as a band and not as a number

    A range against the mean of the bonds the engine did lay out; a literal coordinate would gate the
    bisector formula, which is an implementation.
    """
    if ctx is None:
        skip('quickjs is not installed')

    mol = smiles('F[C@]([H])(Cl)Br')
    plane = mol.layout2d(engine='smilesdrawer')

    hydrogen, = (a.n for a in mol.atoms() if a.atomic_symbol == 'H')
    partner, = (b.m if b.n == hydrogen else b.n for b in mol.bonds()
                if hydrogen in (b.n, b.m))
    heavy = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1])
             for b in mol.bonds() if hydrogen not in (b.n, b.m)]
    mean = sum(heavy) / len(heavy)
    reach = hypot(plane[hydrogen][0] - plane[partner][0], plane[hydrogen][1] - plane[partner][1])
    assert .5 * mean < reach < 1.5 * mean, f'the hydrogen is at {reach / mean:.2f} of a bond length'


def test_both_of_waters_explicit_hydrogens_are_placed():
    """two deferred hydrogens on ONE neighbour, which is where placing them independently would collide"""
    if ctx is None:
        skip('quickjs is not installed')

    mol = smiles('[H]O[H]')
    mol.clean2d(engine='smilesdrawer')

    lengths = _bond_lengths(mol)
    assert len(mol.coordinates()) == 3
    assert _min_separation(mol) / (sum(lengths) / len(lengths)) > .45, 'the two hydrogens collided'


def test_a_molecule_with_no_explicit_hydrogen_lays_out_exactly_as_it_did():
    """the pinned plane of acetic acid, to the last bit

    A REGRESSION PIN and nothing else -- literal coordinates are not a property of the algorithm: these
    four points were measured before explicit hydrogens were deferred out of the parse tree, and the
    deferral must not move a molecule that has none.  When a smilesdrawer upgrade makes this fail,
    re-measure the four points and update the literals here, never widen the tolerance.
    """
    if ctx is None:
        skip('quickjs is not installed')

    plane = smiles('CC(=O)O').layout2d(engine='smilesdrawer')

    assert plane[1] == approx((0., 0.), abs=1e-12)
    assert plane[2] == approx((.824999999998, 2.020263e-06), abs=1e-12)
    assert plane[3] == approx((1.237499999998, -.714468937859), abs=1e-12)
    assert plane[4] == approx((1.237496500795, .714474998639), abs=1e-12)


def test_a_layout_short_of_a_point_is_refused_rather_than_mis_assigned():
    """the guard on `zip(order, xy)`, which turns a silent wrong picture into a named refusal

    `ImplementationError` and not `KeyError`, reported with both counts at the call that can still see
    them rather than several frames later.
    """
    if ctx is None:
        skip('quickjs is not installed')

    real = ctx

    def short(tree, run_gc=False):
        return real(tree, run_gc=run_gc)[:-1]

    mol = smiles('c1ccccc1')
    original = layout_molecule.ctx
    try:
        layout_molecule.ctx = short
        with raises(ImplementationError, match='5 points for 6 atoms'):
            mol.layout2d(engine='smilesdrawer')
    finally:
        layout_molecule.ctx = original
    # and the real engine still answers through the restored name
    assert len(mol.layout2d(engine='smilesdrawer')) == 6


def test_dihydrogen_lays_out():
    """`[H][H]` is a valid two-atom molecule and receives a valid two-point plane

    For a two-atom molecule there is exactly one drawing up to rotation, so both atoms are in the plane,
    not on the same point, and separated by about the rescaled bond length.
    """
    if ctx is None:
        skip('quickjs is not installed')

    mol = smiles('[H][H]')
    plane = mol.layout2d(engine='smilesdrawer')

    assert len(plane) == 2, 'both atoms must be in the plane'
    (x1, y1), (x2, y2) = plane.values()
    sep = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** .5
    assert sep > 0, 'the two atoms must not be on the same point'
    # the rescaled bond length is 0.825; allow a factor of 2 for rounding across rescale
    assert sep == approx(.825, rel=.1), f'separation {sep:.4f} is not a bond length'


def test_hydrogen_only_mixtures_lay_out():
    """every mixture containing a hydrogen-only component is drawable

    The engine gives no index to a hydrogen component chained through `next`, and the component ORDER
    matters, so both orderings are covered: the plane is total and no two atoms share a point.
    """
    if ctx is None:
        skip('quickjs is not installed')

    def _min_sep(plane):
        pts = list(plane.values())
        return min(((pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2) ** .5
                   for i in range(len(pts)) for j in range(i + 1, len(pts)))

    for smi in ('C.[H-]', '[H-].C', 'C.[H][H]', 'O.[H][H]'):
        mol = smiles(smi)
        plane = mol.layout2d(engine='smilesdrawer')
        assert len(plane) == len(mol), f'{smi}: plane is not total ({len(plane)} / {len(mol)} atoms)'
        assert _min_sep(plane) > 0, f'{smi}: two atoms on the same point'


def test_hydrogenation_reaction_depict_draws_hh_bond():
    """`C=C.[H][H]>>CC` returns an SVG and the H-H bond is drawn as one stroke

    Pinned end-to-end against the same reaction minus dihydrogen, which is the only way to count the one
    stroke the H-H bond adds.
    """
    if ctx is None:
        skip('quickjs is not installed')

    import re

    def _stroked(svg):
        return len([tag for tag in re.findall(r'<path[^>]*/>', svg) if 'stroke=' in tag])

    svg_with = smiles('C=C.[H][H]>>CC').depict()
    svg_without = smiles('C=C>>CC').depict()
    assert '<path' in svg_with, 'nothing was drawn'
    assert _stroked(svg_with) == _stroked(svg_without) + 1, (
        'H-H bond must add exactly one stroked path to the reaction picture')


def test_hydride_alone_lays_out_to_origin():
    """`[H-]` as a sole molecule lays out to the same single-point plane it always did

    It now goes through `_place_deferred_hydrogens` rather than the engine as a lone root, and the answer
    has to be unchanged: `{1: (0.0, 0.0)}`, measured at `cf2dbe3`.
    """
    if ctx is None:
        skip('quickjs is not installed')

    plane = smiles('[H-]').layout2d(engine='smilesdrawer')

    assert len(plane) == 1
    (n, (x, y)), = plane.items()
    assert x == approx(0., abs=1e-12)
    assert y == approx(0., abs=1e-12)
