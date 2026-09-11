# -*- coding: utf-8 -*-
"""How many `roles.tsv` rows can cut a stereogenic centre?

The question a `stereo` column would answer is about the CUT ATOM only: a configuration anywhere else
in the fragment is re-based by the arena with nothing stated, so only the atom the cap lands on -- and,
for a cis/trans unit, the double bond it sits on -- can ever need a template to speak about it.

A row is counted capable when some molecule its group matches yields a product whose site carries a
unit `stereogenic_units()` calls stereogenic.  The witness is built by growing the row's own example:
distinct alkyl chains are hung off the site (and off its double-bond partner) until the site is fully
substituted.  Every chain has a different length, so no two substituents are equivalent.  Growth that
breaks the group's own SMARTS makes the row NOT fire, which is the discrimination the count needs --
`alkyl_deoxy/primary_alcohol` stops matching the moment its carbon becomes secondary.

Growth spends an implicit hydrogen per bond and every probe is `check_valence`d, because the count is
worthless if the witness is not a molecule: a plan that would put a third single bond on a carbonyl
oxygen is dropped rather than stored and asked about.
"""
from collections import Counter
from chython import smiles
from chython.chemistry import calc_implicit, check_valence
from chython.core import H_UNKNOWN
from chython.reactions._tables import ROLE_CAP, functional_groups, roles


def plans(product, marker, site, partner):
    """Growth plans for one cut, cheapest first: `[(atom, how many chains)]`.

    Three moves, because a centre is stereogenic for two different reasons.  Growing the SITE raises its
    substitution -- a primary halide becomes a secondary one.  Growing a NEIGHBOUR breaks a symmetry
    without touching the count, which is the only way to reach 2-butanol from isopropanol: the group's
    own SMARTS pins the site's hydrogens, so adding a fourth substituent there stops the row firing.
    Growing the double-bond PARTNER is the cis/trans equivalent -- both ends need two unlike ones.
    """
    yield []
    for k in (1, 2, 3):
        yield [(site, k)]
    neighbours = [n for n in product.neighbors_of(site) if n != marker and n != partner]
    for n in neighbours:
        yield [(n, 1)]
        yield [(site, 1), (n, 1)]
    # Two at once, because a site with no hydrogen to spend can only be reached this way: tert-butanol's
    # three methyls take two different chains before the carbon has three unlike substituents.
    for i, n in enumerate(neighbours):
        for m in neighbours[i + 1:]:
            yield [(n, 1), (m, 1)]
    if partner is not None:
        yield [(partner, 1)]
        yield [(partner, 2)]
        yield [(site, 1), (partner, 1)]
        for n in neighbours:
            yield [(n, 1), (partner, 1)]


def grow(source, plan):
    """A copy of `source` with `count` distinct alkyl chains hung off each atom in `plan`, or None.

    None when an atom has not got `count` hydrogens to spend or the result fails `check_valence`.
    An unknown count spends nothing: what a hydrogen would displace there is not derivable either.

    Chain lengths start ABOVE the seed's atom count, so no added chain can duplicate a substituent the
    seed already had -- a methyl added to a bromoethane's CH2 gives two methyls, and two identical
    substituents are exactly what makes a centre not stereogenic.
    """
    if any(not count <= source.implicit_h_of(atom) < H_UNKNOWN for atom, count in plan):
        return None
    mol = source.copy()
    length = len(source) + 1
    with mol.edit() as e:
        for atom, count in plan:
            for _ in range(count):
                prev = atom
                for _ in range(length):
                    c = e.add_atom('C')
                    e.add_bond(prev, c, 1)
                    prev = c
                length += 1
    # A stored hydrogen count is a statement and `add_bond` does not silently re-derive it, so the atom
    # that spent one says so here.  Only the grown atoms: recomputing the whole probe would overwrite
    # what its own SMILES stated, `[B-](F)(F)F` among it.
    for atom, _ in plan:
        calc_implicit(mol, atom)
    if any(verdict == 'violation' for _, verdict in check_valence(mol)):
        return None
    return mol


def cut(row, mol):
    """Every (product, marker, site, partner) the row yields; partner is a double-bond end or None."""
    for reaction, where in row.template(mol, report=True):
        if ROLE_CAP not in where:
            continue
        marker = where[ROLE_CAP]
        product = next(p for p in reaction.products if marker in p.atom_numbers)
        site = next(iter(product.neighbors_of(marker)))
        partner = next((n for n in product.neighbors_of(site) if product.order_of(site, n) == 2), None)
        yield product, marker, site, partner


def stereogenic_at_site(row, mol):
    """The first product whose site holds a stereogenic unit, or None."""
    for product, _, site, partner in cut(row, mol):
        anchors = {site} if partner is None else {site, partner}
        if any(u['anchor'] in anchors for u in product.stereogenic_units()):
            return product
    return None


def witness(row, example):
    """The smallest grown probe that makes the row's site stereogenic, or None."""
    seed = smiles(example)
    seed.canonicalize()
    shape = next(cut(row, seed), None)
    if shape is None:
        return None, None, 'does not fire on its example'
    for plan in plans(*shape):
        probe = grow(seed, plan)
        if probe is None:
            continue
        probe.canonicalize()
        product = stereogenic_at_site(row, probe)
        if product is not None:
            return probe, product, None
    return None, None, None


KINDS = {0: 'tetrahedral', 1: 'cis/trans', 2: 'allene'}


def unit_kind(row, probe):
    """Which kind of unit the witness's site holds."""
    for product, _, site, partner in cut(row, probe):
        anchors = {site} if partner is None else {site, partner}
        for u in product.stereogenic_units():
            if u['anchor'] in anchors:
                return KINDS.get(u['kind'], u['kind'])
    return '?'


def site_shape(row, example):
    """What the cut atom IS on the row's own example -- the reason a blind row is blind."""
    seed = smiles(example)
    seed.canonicalize()
    shape = next(cut(row, seed), None)
    if shape is None:
        return 'does not fire'
    product, marker, site, _ = shape
    orders = [product.order_of(site, n) for n in product.neighbors_of(site) if n != marker]
    symbol = product.atom(site).atomic_symbol
    if 4 in orders:
        return f'{symbol}, aromatic'
    if 3 in orders:
        return f'{symbol}, triple bond'
    if 2 in orders:
        return f'{symbol}, double bond'
    if symbol != 'C':
        return f'{symbol}, single bonds'
    return f'{symbol}, sp3 with {product.implicit_h_of(site)} H the group pins'


def main():
    known = functional_groups()
    capable, blind, unfired = [], [], []
    for name, rows in roles().items():
        for row in rows:
            example = row.example or known[row.group].example
            probe, product, note = witness(row, example)
            if note:
                unfired.append((row, note))
            elif probe is None:
                blind.append((row, site_shape(row, example)))
            else:
                capable.append((row, probe, product, unit_kind(row, probe)))

    rows_total = sum(len(rows) for rows in roles().values())
    can_roles = {row.name for row, _, _, _ in capable}
    print(f'{rows_total} rows / {len(roles())} roles\n'
          f'  {len(capable)} rows in {len(can_roles)} roles can cut a stereogenic site\n'
          f'  {len(blind)} rows cannot, {len(unfired)} did not fire\n')

    by_kind = Counter(kind for _, _, _, kind in capable)
    for kind, n in by_kind.most_common():
        print(f'  {n:3} rows {kind}')

    print('\nCAN  (role/group, kind, witness -> capped product)')
    for row, probe, product, kind in sorted(capable, key=lambda r: (r[3], r[0].name)):
        print(f'  {kind:11} {row.name}/{row.group:26} {probe}  ->  {product}')

    print('\nCANNOT, by what the cut atom is')
    blind_shapes = Counter(shape for _, shape in blind)
    for shape, n in blind_shapes.most_common():
        names = sorted({row.name for row, s in blind if s == shape})
        print(f'  {n:3} rows  {shape:38} {", ".join(names)}')

    for row, note in unfired:
        print(f'  !! {row.id} {row.name}/{row.group}: {note}')


if __name__ == '__main__':
    main()
