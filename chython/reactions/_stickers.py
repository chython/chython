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
"""Cut a coupling handle off a molecule and cap the cut with an R.

A sticky fragment is a molecule with one attachment point; a sticky linker has two.  Both exist to be
concatenated, and the rule that makes concatenation work is **the joining bond belongs to the LEFT
partner**: the left spelling keeps its leading bond, the right spelling drops its trailing one, so
`A.sticky_right + B.sticky_left` emits exactly one bond token.

The cap is the template's, spelled `[#0:20]` in every `roles.tsv` row, and `report=True` says which atom
the marker landed on.  So the cut centre exists inside the patch and keeps its configuration, and the
site's hydrogen count comes out of the ordinary recompute.  What is left here is the R index, which a
SMIRKS product side has no spelling for.
"""
from typing import Iterator, NamedTuple, Optional
from ._tables import ROLE_CAP, Role, roles


__all__ = ['StickyFragment', 'StickyLinker', 'sticky_fragments', 'sticky_linkers']


class StickyFragment(NamedTuple):
    """One cut, in the three forms a consumer needs.

    `canonical_smiles` is the dedup key: two cuts of two molecules that yield the same fragment yield
    the same string, which is what the R index being part of the canonical record buys.  `atom` is the
    source molecule's atom the cap hangs off, not necessarily the group's `:1`.
    """
    role: str
    sticky_left: str            # `-c1ccccc1`  -- carries the bond token, glues onto a piece before it
    sticky_right: str           # `c(cccc1)c1` -- carries NO bond token; the next piece supplies it
    canonical_smiles: str
    atom: int                   # source atom the cap hangs off; ids survive a patch


class StickyLinker(NamedTuple):
    """One bi-attachment cut.  `canonical_smiles` is always R1 = left, R2 = right.

    `atom_left` and `atom_right` are the source molecule's atoms the two caps hang off, not necessarily
    the groups' `:1`.
    """
    role_left: str
    role_right: str
    sticky_left: str            # `-A...B`, role_left's end first; leading bond token, none trailing
    sticky_right: str           # `-B...A`, the same linker flipped
    canonical_smiles: str
    atom_left: int
    atom_right: int


def _selected(role: Optional[str], present: dict) -> list[Role]:
    """The rows whose handle the molecule actually carries.

    Filtered against `functional_groups()` up front: capping a cut never creates a coupling handle, so
    a handle absent from the source is absent from every intermediate, and running the patch anyway
    only costs time.
    """
    table = roles()
    if role is None:
        items = table.items()
    elif role in table:
        items = [(role, table[role])]
    else:
        raise ValueError(f'unknown role {role!r}; roles.tsv names {len(table)} roles')
    return [row for _, rows in items for row in rows if row.group in present]


def _unmasked(template, molecule, masked: frozenset):
    """`template(molecule, report=True)`, with the mask applied BEFORE the distinct-outcome collapse.

    The collapse keeps the first of two outcomes that build one structure, and cutting either nitrogen
    of `C1CNCCN1` builds one structure; masking the survivor's site would leave nothing.  So a mask
    turns the collapse off and redoes it here, on the same key -- the products' `canonical_bytes` --
    over the outcomes the mask allows.  Unmasked, the template's own collapse is the whole answer.
    """
    if not masked:
        yield from template(molecule, report=True)
        return
    seen = set()
    for reaction, where in template(molecule, report=True, dedupe=False):
        marker = where[ROLE_CAP]
        # The product HOLDING the marker, not the first one: a row whose patch releases its leaving
        # group as a molecule rather than deleting it yields two products in an unspecified order.
        product = next(p for p in reaction.products if marker in p.atom_numbers)
        if next(iter(product.neighbors_of(marker))) in masked or not masked <= set(product):
            continue
        key = tuple(sorted(p.canonical_bytes for p in reaction.products))
        if key not in seen:
            seen.add(key)
            yield reaction, where


def _index(product, marker, index: int) -> None:
    """Number one end of a linker.  A fragment's single marker keeps index 0 and needs no write."""
    with product.edit() as e:
        e.set_r_index(marker, index)


def sticky_fragments(molecule, role=None, *, masked=None,
                     hydrogens: bool = False) -> Iterator[StickyFragment]:
    """Every mono-attachment fragment this molecule exposes, one per match.

    `masked` bars an atom from the coupling in both the ways it can take part: as the attachment site,
    and as a leaving group the patch consumes.  Atom ids survive a patch, so a masked atom is consumed
    exactly when it is absent from the product.  An equivalent unmasked atom still yields:
    `C1CNCCN1` masked at one nitrogen is cut at the other.

    A molecule with more than one connected component yields nothing -- a counter-ion has no cut.  A
    salt is therefore the caller's to reduce to one component: `decompose_salts().parents` names them, one
    row per component, and the anion of a trifluoroborate or a carboxylate enumerates normally.

    A fragment's marker has no valence rules, so `check_valence` reports it as unknown; a caller
    filtering on a clean valence report excludes the marker.
    """
    if molecule.connected_components_count != 1:
        return
    masked = frozenset(masked or ())
    for row in _selected(role, molecule.functional_groups()):
        for reaction, where in _unmasked(row.template, molecule, masked):
            marker = where[ROLE_CAP]
            product = next(p for p in reaction.products if marker in p.atom_numbers)
            site = next(iter(product.neighbors_of(marker)))
            product.canonicalize()
            yield StickyFragment(
                row.name,
                product.sticky_smiles(left=marker, remove_left=True, keep_bond_left=True,
                                      hydrogens=hydrogens),
                product.sticky_smiles(right=marker, remove_right=True, hydrogens=hydrogens),
                str(product), site)


def sticky_linkers(molecule, role_left=None, role_right=None, *, masked=None,
                   hydrogens: bool = False) -> Iterator[StickyLinker]:
    """Every bi-attachment linker this molecule exposes, one per (left match, right match).

    `masked` applies to the LEFT end only, by design: a masked handle is one whose only role is the
    deferred second step, so it may sit on the right and never on the left.  The asymmetry is also what
    removes the (masked_left, X) / (X, masked_right) ordering duplicate.

    A linker whose two caps would land on the same atom is skipped: a linker needs at least one atom
    between its ends.

    A molecule with more than one connected component yields nothing -- a counter-ion has no cut.  A
    salt is therefore the caller's to reduce to one component: `decompose_salts().parents` names them, one
    row per component, and the anion of a trifluoroborate or a carboxylate enumerates normally.
    """
    if molecule.connected_components_count != 1:
        return
    masked = frozenset(masked or ())
    present = molecule.functional_groups()
    left_rows = _selected(role_left, present)
    right_rows = _selected(role_right, present)
    if not left_rows or not right_rows:
        return

    for left_row in left_rows:
        for left_reaction, left_where in _unmasked(left_row.template, molecule, masked):
            one = left_where[ROLE_CAP]
            inter = next(p for p in left_reaction.products if one in p.atom_numbers)
            # read on `inter`: the second cut cannot rebond this site, since `one` must survive it
            site_one = next(iter(inter.neighbors_of(one)))
            _index(inter, one, 1)
            for right_row in right_rows:
                for reaction, where in right_row.template(inter, report=True):
                    two = where[ROLE_CAP]
                    product = next(p for p in reaction.products if two in p.atom_numbers)
                    # The first end has to survive the second cut, on the same fragment: a right row
                    # whose leaving group swallows it yields a mono-attachment fragment, not a linker.
                    if one not in product.atom_numbers:
                        continue
                    # Both caps on one atom: a linker needs at least one atom between its ends.
                    site_two = next(iter(product.neighbors_of(two)))
                    if site_two == site_one:
                        continue
                    _index(product, two, 2)
                    product.canonicalize()
                    yield StickyLinker(
                        left_row.name, right_row.name,
                        product.sticky_smiles(left=one, right=two, remove_left=True,
                                              keep_bond_left=True, remove_right=True,
                                              hydrogens=hydrogens),
                        product.sticky_smiles(left=two, right=one, remove_left=True,
                                              keep_bond_left=True, remove_right=True,
                                              hydrogens=hydrogens),
                        str(product), site_one, site_two)
