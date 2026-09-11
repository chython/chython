# -*- coding: utf-8 -*-
#
#  Copyright 2019-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019, 2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
"""2D layout of one molecule, as module functions over the core container.

Every entry point takes `mol` first; `chython.depict._hooks` injects the public ones onto the container.
A plane is a plain `{n: (x, y)}` dict and is never the molecule -- only `clean2d`, `rescale2d`
and `_store_plane` write to the arena, so a renderer can lay a molecule out without changing it.
"""
from importlib.resources import files
from math import atan2, cos, fsum, hypot, pi, sin
from ...exceptions import ImplementationError
from .._config import Clean2DEngine, get_clean2d_engine

# `clean2d.js` is an esbuild IIFE bundle publishing a single global `$`; this line gives the binding a
# top-level name to fetch, and lives here so the shipped bundle stays byte-identical to esbuild's output.
_JS_SHIM = '\nfunction clean2d(tree) { return $.clean2d(tree); }\n'

# smiles-drawer's own `bondLength`, in the units its layout comes back in.  Used as the spacing for a
# deferred hydrogen when the engine's plane holds no bond length at all (`C[H]`, `[H]O[H]`), so that
# such a molecule keeps one scale before `_rescale_plane` normalizes it.
_ENGINE_BOND_LENGTH = 30.

try:
    from quickjs import Function as _JSFunction

    # `Function`, not `Context`: it pins the QuickJS runtime to one dedicated worker thread behind a
    # lock, where a bare `Context` crashes when touched from a second thread even sequentially -- and
    # this singleton is shared by every caller.  It also marshals the parse tree through QuickJS's JSON.
    ctx = _JSFunction('clean2d', files(__package__).joinpath('clean2d.js').read_text(encoding='utf-8') + _JS_SHIM)
except Exception:  # absent, or built against an incompatible libquickjs
    ctx = None


def _adjacency(mol):
    """`{n: {neighbour: order}}` from one pass over `mol.bonds()`, built once because the walks
    below query every bond several times and `order_of` is a binary search.

    Every atom gets a row, isolated ones included.  Rows are in arena order and each row in `bonds()`
    order: a set anywhere here would make the spanning forest depend on hash order, and the layout too.
    """
    out = {}
    for n in mol:
        out[n] = {}
    for bond in mol.bonds():
        out[bond.n][bond.m] = bond.order
        out[bond.m][bond.n] = bond.order
    return out


def _stored_plane(mol):
    """The molecule's own coordinates as a plane, with the origin standing in for "none stated".

    `coordinates()` answers `{}` when the arena carries no XY segment; the geometry below subscripts the
    plane for every atom of every bond, so it needs a map total over `mol`.
    """
    plane = mol.coordinates()
    if plane:
        return plane
    return dict.fromkeys(mol, (0., 0.))


def layout2d(mol, *, engine: Clean2DEngine = None, force: bool = False):
    """Compute a 2d layout and return it as `{n: (x, y)}`, leaving the molecule untouched.

    This is the form a renderer wants -- drawing must not change what it draws.  `clean2d()` is this
    plus the decision to keep the result.  `force=False` on a molecule that already `has_layout`
    returns the stored plane and computes nothing.  By default the JS implementation of
    https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425 is used; it can be changed globally with the
    `chython.clean2d_engine` parameter.

    :param engine: override globally set engine
    :param force: recompute even if the molecule already carries a layout
    """
    if not force and mol.has_layout:
        return mol.coordinates()

    plane = _engine_layout(mol, get_clean2d_engine(engine))
    _rescale_plane(mol, plane)
    if mol.connected_components_count > 1:
        shift_x = 0.
        for c in mol.connected_components:
            shift_x = _shift_plane_mean(mol, plane, shift_x, component=c) + .9
    return plane


def clean2d(mol, *, engine: Clean2DEngine = None, force: bool = False):
    """Compute a 2d layout and store it on the molecule.

    Not always a recomputation: a molecule that already `has_layout` is left exactly as it is, since the
    request is "make sure this molecule has a layout".  `force=True` relays it regardless.

    :param engine: override globally set engine
    :param force: lay the molecule out again whatever coordinates it already has
    """
    if not force and mol.has_layout:
        return
    # `force=True` below: the has_layout question is already answered, and asking it again would send a
    # molecule that has a layout down the stored-plane branch.
    _store_plane(mol, layout2d(mol, engine=engine, force=True))


def _engine_layout(mol, engine: Clean2DEngine):
    """`{n: (x, y)}` from the named backend, unrescaled and unshifted.

    Every `to_*` import sits inside its branch: an unnamed toolkit must not be imported to lay a
    molecule out, and `chython.depict` may not import `chython.interop` at module scope.
    """
    plane = {}
    if engine == 'rdkit':
        from rdkit.Chem.AllChem import Compute2DCoords
        from ...interop._rdkit import to_rdkit

        rd = to_rdkit(mol, keep_mapping=False)
        Compute2DCoords(rd)
        # set coordinates from the first rdkit conformer. usually it's 2d layout
        for n, (x, y, _) in zip(mol, rd.GetConformers()[0].GetPositions()):
            plane[n] = (x, y)
    elif engine == 'smilesdrawer':
        if ctx is None:
            raise ImportError('quickjs is not installed or broken')
        # smiles-drawer normalizes the layout regardless of the tree root, so a single
        # deterministic layout pass is enough.
        tree, order = _clean2d_tree(mol)
        if not order:
            # Every atom was withheld by `_deferred_hydrogens`: no component has a heavy atom
            # (`[H][H]`, `[H-]`).  The engine is not called -- there is no tree and `xy` would be empty.
            _place_deferred_hydrogens(mol, plane)
            return plane
        try:
            # `run_gc=False`: the binding otherwise runs a full QuickJS collection after every call,
            # which costs 3x and reclaims nothing -- the engine's own threshold GC holds this context
            # under a megabyte across thousands of layouts.
            xy = ctx(tree, run_gc=False)
        except Exception:
            raise ImplementationError

        # The `zip` below is a positional correspondence, and `zip` truncates in silence: with fewer
        # points than atoms every atom after the missing one takes its neighbour's place, so the
        # molecule is laid out wrong rather than not laid out.  Check it here instead.
        if len(xy) != len(order):
            raise ImplementationError(f'smiles-drawer returned {len(xy)} points for {len(order)} '
                                      f'atoms: the layout cannot be assigned')
        # `float()` and not the bare subtraction: QuickJS returns an integral JS number as a Python
        # `int`, so a one-atom molecule laid out at the origin would come back `(0, 0)`.
        shift_x, shift_y = xy[0]
        for n, (x, y) in zip(order, xy):
            plane[n] = (float(x - shift_x), float(shift_y - y))
        # The explicit hydrogens the tree left out, here and not later: everything downstream
        # subscripts the plane per atom, so the plane this branch returns must be total over `mol`.
        _place_deferred_hydrogens(mol, plane)
    elif engine == 'cdk':
        from ...interop._cdk import to_cdk
        from ...interop._java import get_cdk

        sdg = get_cdk().layout.StructureDiagramGenerator()
        sdg.setUseTemplates(False)
        sdg.setMolecule(to_cdk(mol))
        sdg.generateCoordinates()
        cdk_mol = sdg.getMolecule()

        # to_cdk preserves atom order: CDK atom index i (0-based) matches the i-th atom
        for i, n in enumerate(mol):
            xy = cdk_mol.getAtom(i).getPoint2d()
            plane[n] = (xy.x, xy.y)
    elif engine == 'obabel':
        from openbabel import openbabel
        from ...interop._openbabel import to_openbabel

        ob = to_openbabel(mol)
        assert openbabel.OBOp.FindType('gen2D').Do(ob), 'OpenBabel failed to generate 2d layout'
        assert ob.NumAtoms() == len(mol), 'OpenBabel modified molecule'

        # to_openbabel preserves atom order: OBMol index i (1-based) matches the i-th atom
        for i, n in enumerate(mol, 1):
            xy = ob.GetAtom(i).GetVector()
            plane[n] = (xy.GetX(), xy.GetY())
    elif engine == 'indigo':
        from ...interop._indigo import to_indigo

        ind = to_indigo(mol)
        assert not ind.layout(), 'Indigo failed to generate 2d layout'

        # to_indigo preserves atom order: iterateAtoms() matches mol order
        for n, a in zip(mol, ind.iterateAtoms()):
            x, y, _ = a.xyz()
            plane[n] = (x, y)
    else:
        raise ValueError(f'Invalid clean2d engine: {engine}')
    return plane


# The geometry below operates on a plane, never on the molecule; `rescale2d` and the `_fix_plane_*`
# functions are the same operations applied to the stored coordinates.


def _store_plane(mol, plane):
    """Write a plane into the arena, in one edit scope -- leaving the scope rebuilds it, so a scope per
    atom would rebuild it once per coordinate."""
    with mol.edit():
        for n, (x, y) in plane.items():
            mol.set_xy(n, x, y)


def _rescale_plane(mol, plane) -> bool:
    bonds = []
    for bond in mol.bonds():
        nx, ny = plane[bond.n]
        mx, my = plane[bond.m]
        bonds.append(hypot(nx - mx, ny - my))
    if bonds:
        bond_reduce = fsum(bonds) / len(bonds) / .825
        if bond_reduce > .5:  # check for singularity
            for n, (x, y) in plane.items():
                plane[n] = (x / bond_reduce, y / bond_reduce)
            return True
    return False


def _shift_plane_mean(mol, plane, shift_x: float, shift_y=0., component=None) -> float:
    if component is None:
        component = plane

    left = min(component, key=lambda x: plane[x][0])
    right = max(component, key=lambda x: plane[x][0])

    min_x = plane[left][0] - shift_x
    if len(mol.atom(left).atomic_symbol) == 2:
        min_x -= .2

    max_x = plane[right][0] - min_x
    min_y = min(plane[x][1] for x in component)
    max_y = max(plane[x][1] for x in component)
    mean_y = (max_y + min_y) / 2 - shift_y
    for n in component:
        x, y = plane[n]
        plane[n] = (x - min_x, y - mean_y)

    if -.18 <= plane[right][1] <= .18:
        # `implicit_h` is None for an undeterminable count, which is falsy and so takes the "no
        # hydrogens" branch: an unknown count must not widen the box by a label that may not be drawn.
        factor = mol.atom(right).implicit_h
        if factor == 1:
            max_x += .15
        elif factor:
            max_x += .25
    return max_x


def _shift_plane_min(mol, plane, shift_x: float, shift_y=0., component=None) -> float:
    if component is None:
        component = plane

    right = max(component, key=lambda x: plane[x][0])
    min_x = min(plane[x][0] for x in component) - shift_x
    max_x = plane[right][0] - min_x
    min_y = min(plane[x][1] for x in component) - shift_y
    for n in component:
        x, y = plane[n]
        plane[n] = (x - min_x, y - min_y)

    if shift_y - .18 <= plane[right][1] <= shift_y + .18:
        factor = mol.atom(right).implicit_h
        if factor == 1:
            max_x += .15
        elif factor:
            max_x += .25
    return max_x


def rescale2d(mol) -> bool:
    """Rescale stored coordinates to average bond length 0.825, and answer whether it rescaled.

    False, and nothing stored, when there is no scale to read: no coordinates, no bonds, or a plane
    collapsed tightly enough that dividing by its mean would be a singularity.  Registered onto
    `MoleculeContainer.rescale2d`, whose docstring is the one a caller reads.
    """
    plane = _stored_plane(mol)
    if _rescale_plane(mol, plane):
        _store_plane(mol, plane)
        return True
    return False


def _fix_plane_mean(mol, shift_x: float, shift_y=0., component=None) -> float:
    plane = _stored_plane(mol)
    max_x = _shift_plane_mean(mol, plane, shift_x, shift_y, component)
    _store_plane(mol, plane)
    return max_x


def _fix_plane_min(mol, shift_x: float, shift_y=0., component=None) -> float:
    plane = _stored_plane(mol)
    max_x = _shift_plane_min(mol, plane, shift_x, shift_y, component)
    _store_plane(mol, plane)
    return max_x


def _deferred_hydrogens(mol, bonds):
    """The atoms left out of the parse tree, to be placed once the engine has answered.

    smiles-drawer gives an atom an index only when it is not H or is a lone root, so an explicit `[H]`
    in the tree costs a point that never comes back.  Read by `_place_deferred_hydrogens` too, so the
    two cannot drift.  Withheld: every atom of a component with no heavy atom (`[H]`, `[H-]`, `[H][H]`,
    with no lone-root exemption -- the engine's own fails at the first chained component, as `C.[H-]`
    shows), and a hydrogen of degree exactly one whose neighbour is heavy.  Degree two or more stays in
    the tree, since a bridging hydride would take a spanning-forest edge with it.
    """
    out = []
    # Group 1: every atom in a hydrogen-only component.
    deferred_set = set()
    for component in mol.connected_components:
        if all(mol.atom(n).atomic_symbol == 'H' for n in component):
            out.extend(component)
            deferred_set.update(component)
    # Group 2: H of degree 1 beside a heavy neighbour.
    for atom in mol.atoms():
        sid = atom.n
        if sid in deferred_set or atom.atomic_symbol != 'H' or len(bonds[sid]) != 1:
            continue
        partner, = bonds[sid]
        if mol.atom(partner).atomic_symbol != 'H':
            out.append(sid)
    return out


def _place_deferred_hydrogens(mol, plane):
    """Put every atom `_deferred_hydrogens` withheld back into the plane.

    1. a component with no heavy atom: a chain along +x from the lowest stable id.  The origin is a fine
       start because `layout2d` shifts the components apart afterwards.
    2. an explicit hydrogen beside a heavy neighbour: the widest angular gap between that neighbour's
       other bonds.  One at a time in stable-id order, so a second hydrogen on one atom sees the first
       as occupied -- which is what places water's two.

    Spacing is the mean engine-placed bond length, so `_rescale_plane` normalises at one scale.
    """
    bonds = _adjacency(mol)
    lengths = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1])
               for b in mol.bonds() if b.n in plane and b.m in plane]
    fallback = fsum(lengths) / len(lengths) if lengths else _ENGINE_BOND_LENGTH

    # Case 1: hydrogen-only components.
    for component in mol.connected_components:
        if all(mol.atom(n).atomic_symbol == 'H' for n in component):
            for i, n in enumerate(sorted(component)):
                plane[n] = (i * fallback, 0.)

    # Case 2: H atoms beside heavy neighbours.
    for n in sorted(_deferred_hydrogens(mol, bonds)):
        if n in plane:
            continue  # already placed in case 1
        partner, = bonds[n]
        px, py = plane[partner]
        taken = [(atan2(plane[other][1] - py, plane[other][0] - px),
                  hypot(plane[other][0] - px, plane[other][1] - py))
                 for other in bonds[partner] if other in plane]
        if taken:
            angles = sorted(angle for angle, _ in taken)
            distance = fsum(length for _, length in taken) / len(taken)
            # the widest circular gap, closing round through the first angle.  One neighbour gives one
            # gap of 2*pi, hence the direction straight opposite it.
            width, start = max((angles[(i + 1) % len(angles)] - angle + (2. * pi if i + 1 == len(angles)
                                                                         else 0.), angle)
                               for i, angle in enumerate(angles))
            direction = start + width / 2.
        else:
            distance, direction = fallback, 0.
        plane[n] = (px + distance * cos(direction), py + distance * sin(direction))


def _clean2d_tree(mol):
    """Build a smiles-drawer parse tree from the molecule graph -- a plain iterative DFS with
    ring-closure detection, since the layout needs only valid connectivity and is root-invariant.

    Returns `(tree, order)` where `order[i]` is the stable id of the i-th atom created, matching
    smiles-drawer's atom index; the `_deferred_hydrogens` are in neither.  See clean2d/README for the
    node schema.
    """
    atoms = {a.n: a for a in mol.atoms()}
    bonds = _adjacency(mol)
    # Filtering the adjacency once means the DFS, the branch lists and the ring-closure pass all stop
    # seeing the deferred hydrogens at once.  `bonds` keeps `_adjacency`'s arena order.
    deferred = set(_deferred_hydrogens(mol, bonds))
    bonds = {n: {m: order for m, order in row.items() if m not in deferred}
             for n, row in bonds.items() if n not in deferred}
    order = []
    nodes = {}  # atom number -> its node dict
    parent_of = {}  # atom number -> tree parent atom number (root: None)

    # Layout depends only on element and connectivity, so a node's `atom` is a bare element string.
    # Orders with no single/double/triple counterpart (aromatic 4, dative 8) collapse to '-'.
    bond_symbol = {1: '-', 2: '=', 3: '#', 4: '-', 8: '-'}

    # spanning forest via iterative DFS. all neighbours become `branches`.
    components = []
    for root in bonds:
        if root in nodes:
            continue
        nodes[root] = rnode = {'atom': atoms[root].atomic_symbol, 'isBracket': False,
                               'branches': [], 'branchCount': 0, 'ringbonds': [], 'ringbondCount': 0,
                               'bond': '-', 'branchBond': '-', 'next': None, 'hasNext': False}
        order.append(root)
        components.append(rnode)
        parent_of[root] = None
        stack = [root]
        while stack:
            parent = stack[-1]
            for child in bonds[parent]:
                if child not in nodes:
                    bs = bond_symbol[bonds[parent][child]]
                    nodes[child] = cnode = {'atom': atoms[child].atomic_symbol, 'isBracket': False,
                                            'branches': [], 'branchCount': 0, 'ringbonds': [],
                                            'ringbondCount': 0, 'bond': bs, 'branchBond': bs,
                                            'next': None, 'hasNext': False}
                    order.append(child)
                    pnode = nodes[parent]
                    pnode['branches'].append(cnode)
                    pnode['branchCount'] += 1
                    parent_of[child] = parent
                    stack.append(child)
                    break
            else:
                stack.pop()

    # ring closures: every non-tree edge gets a matching ringbond id on both ends.
    cycle = 0
    for n in order:
        for m in bonds[n]:
            if m <= n or parent_of.get(m) == n or parent_of.get(n) == m:
                continue
            cycle += 1
            bs = bond_symbol[bonds[n][m]]
            for k in (n, m):
                nodes[k]['ringbonds'].append({'bond': bs, 'id': cycle})
                nodes[k]['ringbondCount'] += 1

    if not components:
        # No heavy atoms at all; `_engine_layout` checks `not order` and skips the engine.
        return None, []

    # chain disconnected components through `next` with a '.' bond.
    for prev, comp in zip(components, components[1:]):
        comp['bond'] = '.'
        prev['next'] = comp
        prev['hasNext'] = True
    return components[0], order


__all__ = ['clean2d', 'layout2d', 'rescale2d']
