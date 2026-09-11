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
"""MACCS structural keys, Durant, Leland, Henry, Nourse, J. Chem. Inf. Comput. Sci. 2002, 42, 1273.

The 166-bit numbering is MDL's; the patterns are chython's reading of the published key descriptions,
in `tables/maccs.tsv`.  Nothing here claims bit-for-bit parity with any other implementation, and a bit
that differs from another toolkit's is a documented difference -- widely used implementations knowingly
differ from the published list, and where they do this table follows the publication.  The oracle is
`tables/maccs_corpus.tsv`, whose expected answers are read off the published descriptions by hand.
"""
from ._tables import MACCS_PREDICATES, maccs_rules
from ..core._core import require_numpy


def _has_isotope(molecule) -> bool:
    """Key 1.  Any atom whose isotope is stated -- an atom field, not a substructure."""
    return any(a.isotope for a in molecule.atoms())


def _atomic_number_gt_103(molecule) -> bool:
    """Key 2.  Any atom past lawrencium.

    `atom.element` IS the atomic number; there is no `atom.atomic_number` in chython 3.  Written as a
    comparison rather than a fifteen-element `,` list, because a range is what the key says.
    """
    return any(a.element > 103 for a in molecule.atoms())


def _has_charge(molecule) -> bool:
    """Key 49.  Any atom with a non-zero formal charge."""
    return any(a.charge for a in molecule.atoms())


def _fragments_gt_1(molecule) -> bool:
    """Key 166.  More than one connected component -- a property of the RECORD, not of a fragment."""
    return len(molecule.connected_components) > 1


def _ring_present(molecule) -> bool:
    """Key 165.  At least one ring."""
    return bool(molecule.sssr)


def aromatic_ring_count(molecule) -> int:
    """How many aromatic rings the molecule has.

    A DELEGATION, not a second answer: `MoleculeContainer.aromatic_rings_count` is computed in
    `core/_descriptors.pxi`, and recomputing it from `sssr` plus per-atom hybridization would put a
    second aromatic-ring answer in the tree.  The wrapper exists so that key 125's predicate and QED's
    AROM term name one function.
    """
    return molecule.aromatic_rings_count


def _aromatic_rings_gt_1(molecule) -> bool:
    """Key 125, `Aromatic Ring > 1`.  RINGS, not aromatic atoms -- benzene has six of the latter."""
    return aromatic_ring_count(molecule) > 1


def _six_rings_gt_1(molecule) -> bool:
    """Key 145, `6M Ring > 1`.  Counted over `sssr`, because `[!#1;*;r6]` counts ATOMS: it matches
    benzene six times, so a `count = 2` row would set the key on a single ring."""
    return sum(1 for ring in molecule.sssr if len(ring) == 6) > 1


#: The non-substructure keys, name -> `(molecule) -> bool`.
MACCS_PREDICATE_FNS = {'isotope': _has_isotope,
                       'atomic_number_gt_103': _atomic_number_gt_103,
                       'charge': _has_charge,
                       'fragments_gt_1': _fragments_gt_1,
                       'ring_present': _ring_present,
                       'aromatic_rings_gt_1': _aromatic_rings_gt_1,
                       'six_rings_gt_1': _six_rings_gt_1}

# asserted at IMPORT, not in a test: a name added to the table's vocabulary without a function here
# would otherwise be a load-time `KeyError` on whichever molecule first reached that row.
assert set(MACCS_PREDICATE_FNS) == set(MACCS_PREDICATES), \
    'every registered MACCS predicate name needs a function and vice versa'


def _distinct_matches(query, molecule) -> int:
    """How many DISTINCT atom sets the query matches.

    A count key asks "how many of these are there", and a symmetric pattern maps onto one site several
    ways -- the same double count a rotatable-bond pass gets wrong by counting mappings.  The atom SET
    is the site.
    """
    return len({frozenset(mapping.values()) for mapping in query.get_mapping(molecule)})


def maccs_match_counts(molecule) -> dict[int, int]:
    """Key number -> number of distinct matched atom sets, for the keys that matched at all.

    Predicate keys report 1 when true and are absent when false.  For diagnosing a `maccs_corpus.tsv`
    failure without re-deriving the pattern by hand.
    """
    out = {}
    for row in maccs_rules():
        if row.kind == 'unset':
            continue                      # no definition to match; the bit is permanently zero
        elif row.kind == 'predicate':
            if MACCS_PREDICATE_FNS[row.predicate](molecule):
                out[row.key] = 1
        else:
            n = _distinct_matches(row.query, molecule)
            if n:
                out[row.key] = n
    return out


def maccs_keys(molecule):
    """The 166 published MACCS structural keys as `uint8[167]`.

    ONE-BASED: `keys[n]` is published key `n` for `n` in 1..166, and index 0 is permanently zero so
    that no caller writes `n - 1`.  Key 44 is permanently zero too, and its row says why.

    NUMPY IS IMPORTED HERE AND NOT AT MODULE LEVEL.  `chython.chemistry.__init__` imports this module,
    so a module-level `from numpy import zeros` makes numpy a hard dependency of the whole base install
    -- and `require_numpy()` runs first so the failure is the core's one message naming the `ml` extra
    rather than a bare "No module named 'numpy'" from the line below it.
    """
    require_numpy()
    from numpy import uint8, zeros

    out = zeros(167, dtype=uint8)
    for row in maccs_rules():
        if row.kind == 'unset':
            continue                      # see MACCS_UNSET_KEYS: no stated definition, so no bit
        elif row.kind == 'predicate':
            if MACCS_PREDICATE_FNS[row.predicate](molecule):
                out[row.key] = 1
        elif _distinct_matches(row.query, molecule) >= row.count:
            out[row.key] = 1
    return out


def maccs_bit_set(molecule) -> frozenset:
    """The set of published MACCS key numbers this molecule sets.  1..166, never 0."""
    v = maccs_keys(molecule)
    return frozenset(n for n in range(1, 167) if v[n])
