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
"""Turn a :class:`~chython.formats.pdb._records.PDBRecord` into a container by lookup, never by
perception: the file states a component id and `tables/residues.tsv` says what that component is.

Nothing calls this -- the readers reach for no template -- and one record is one model.  Two choices: a
zero charge yields to the table (a record cannot tell a stated zero from a blank column) while a non-zero
file charge wins, and consecutive numbering is read as the file's statement of chain order.  A finding
that can repeat is counted by kind and reported once with an example; `unsupported: ` means our table is
the limitation and any other prefix means the file was broken and we read it anyway.
"""

from collections.abc import Sequence

from ._records import PDBRecord
from ...core import LogRecord, LOST
from ...core._core import CONF_EXT_INDEX_MAX, MoleculeContainer


__all__ = ['build_molecule']


class _Findings:
    """Findings counted by kind rather than logged residue by residue, keeping the first residue seen.

    The key is `(prefix, what)`, so two phrasings never share a count.  Report order is first-seen,
    which for one record is file order and therefore deterministic.
    """
    __slots__ = ('counts', 'first')

    def __init__(self):
        self.counts: dict[tuple[str, str], int] = {}
        self.first: dict[tuple[str, str], str] = {}

    def hit(self, prefix: str, what: str, where: str) -> None:
        key = (prefix, what)
        self.counts[key] = self.counts.get(key, 0) + 1
        self.first.setdefault(key, where)

    def report(self, log: list[str]) -> None:
        for (prefix, what), count in self.counts.items():
            log.append(LogRecord('pdb:findings', (),
                                 f'{prefix}: {count} {what} (first {self.first[(prefix, what)]})'))


def _label(key: tuple) -> str:
    """The residue a finding points at, as `chain/NAME seq`.  One record is one model, so the model
    number is constant and left out."""
    _, chain, name, sequence, ins = key
    return f'{chain or "?"}/{name or "?"}{"" if sequence is None else sequence}{ins or ""}'


def _match_key(atom) -> tuple:
    """What makes two records' atoms the same atom across models.

    `PDBAtom.residue_key` with its model number dropped, plus the atom's own name and alternate id.
    `serial` is not part of it: no spec pins a serial across models.
    """
    _, chain, name, sequence, ins = atom.residue_key
    return chain, sequence, ins, name, atom.atom_name, atom.alt_loc


def _model_label(record: PDBRecord) -> str:
    return f'model {record.model}' if record.model is not None else 'a model with no number'


def _ext_index(record: PDBRecord, findings: _Findings) -> int | None:
    """This record's model number as an `ext_index`, or None when the field cannot hold it.

    A `MODEL` serial is a signed integer as far as the reader is concerned, so a negative one is
    storable as a model and not as a number; the geometry still lands.
    """
    number = record.model
    if number is None or 0 <= number <= CONF_EXT_INDEX_MAX:
        return number
    findings.hit('conformer', 'model(s) state a number outside the range the conformer record holds, '
                              'so the model is stored and states no number', _model_label(record))
    return None


def _residue_groups(record: PDBRecord) -> dict[tuple, list[int]]:
    """`residue_key` -> the record indexes of its atoms, in file order.

    `PDBAtom.residue_key` leaves `alt_loc` out, so alternate conformers group together here and are
    separated only when one of them is chosen.
    """
    groups: dict[tuple, list[int]] = {}
    for index, atom in enumerate(record.atoms):
        groups.setdefault(atom.residue_key, []).append(index)
    return groups


def _select_conformers(record: PDBRecord, groups: dict[tuple, list[int]], requested: str | None,
                       findings: _Findings) -> set[int]:
    """The record indexes that survive alternate-conformer selection.

    A molecular graph cannot hold one atom at two positions, so one conformer is chosen: every atom with
    no `alt_loc` is kept, plus exactly one alternate id.  The default is the highest summed occupancy,
    ties broken by the smallest id -- not file order, which is a formatting convention and would give
    two molecules for two writers' orderings of one structure.  An unstated occupancy adds nothing.
    """
    kept: set[int] = set()
    for key, indexes in groups.items():
        by_alt: dict[str | None, list[int]] = {}
        for index in indexes:
            by_alt.setdefault(record.atoms[index].alt_loc, []).append(index)
        alternates = sorted(alt for alt in by_alt if alt is not None)
        if not alternates:
            kept.update(indexes)
            continue

        kept.update(by_alt.get(None, ()))
        if requested is None:
            chosen = min(alternates,
                         key=lambda alt: (-sum(record.atoms[i].occupancy or 0.0
                                               for i in by_alt[alt]), alt))
            findings.hit('residue', f'residue(s) hold alternate conformers; the one with the highest '
                                    f'summed occupancy is kept and the rest are dropped',
                         f'{_label(key)} keeps {chosen!r} of ' + ', '.join(repr(a) for a in alternates))
            kept.update(by_alt[chosen])
        elif requested in alternates:
            kept.update(by_alt[requested])
        else:
            findings.hit('residue', f'residue(s) hold alternate conformers but not the requested '
                                    f'{requested!r}, so only their conformer-free atoms are kept',
                         f'{_label(key)} holds ' + ', '.join(repr(a) for a in alternates))
    return kept


class _Residue:
    """One residue that matched a template, as the plan pass left it.

    `by_name` is normalised template atom name -> record index, through which the template bonds and the
    chain link resolve.  `heavy` is every non-hydrogen atom kept, template-named or not.  `stated_h` is
    how many explicit hydrogens this pass dropped.  `missing` is computed by the plan pass and reported
    only after the chain link pass, which is when `links_built` -- the link atom names actually used --
    is known and the shortfall can be suppressed correctly.
    """
    __slots__ = ('key', 'template', 'by_name', 'heavy', 'stated_h', 'missing', 'links_built')

    def __init__(self, key: tuple, template):
        self.key = key
        self.template = template
        self.by_name: dict[str, int] = {}
        self.heavy: list[int] = []
        self.stated_h = 0
        self.missing: list[str] = []
        self.links_built: set[str] = set()


def _plan_templated(record: PDBRecord, residue: _Residue, indexes: list[int],
                    plan: dict[int, tuple], dropped_h: set[int], findings: _Findings) -> None:
    """Decide what to store for one residue the table has a row for.

    Explicit hydrogens are dropped here because `ResidueTemplate.atoms` is heavy atoms only, so one
    would have no template bond and arrive isolated; the count is kept and checked against the
    derivation afterwards.  An atom the template does not name is kept and simply gets no template bond.
    """
    from ...chemistry import normalize_atom_name

    label = _label(residue.key)
    name = residue.key[2]
    kind = residue.template.kind
    for index in indexes:
        atom = record.atoms[index]
        if atom.element == 'H':
            # Element and never atom name: a haem's `NA` is a pyrrole nitrogen, so the name column is
            # not evidence about an element in this format.
            residue.stated_h += 1
            dropped_h.add(index)
            continue

        if atom.atom_name is None:
            row = None
        else:
            # `kind` is passed because the alias map is nucleotide-only: `O1P` is a legacy spelling of
            # `OP1` in a nucleotide and a live atom name in a phosphorylated residue, so normalising
            # before the row is known would rename a real atom.
            key_name = normalize_atom_name(atom.atom_name, kind)
            row = residue.template.atoms.get(key_name)

        if row is None:
            if atom.atom_name is None:
                findings.hit('atom', f'atom(s) of a residue matching a template state no atom name, '
                                     f'so they cannot be joined to it and are stored unbonded', label)
            else:
                findings.hit('residue', f'atom(s) of {name} carry a name the {name} template does not '
                                        f'have ({key_name}); they are stored and get no template '
                                        f'bond', label)
            spec = _atom_spec(atom, None, label, findings)
            if spec is not None:
                plan[index] = spec
                residue.heavy.append(index)
            continue

        if key_name in residue.by_name:
            findings.hit('residue', f'atom(s) of {name} restate the atom name {key_name} within one '
                                    f'conformer; the first is joined to the template and the rest are '
                                    f'stored unbonded', label)
            spec = _atom_spec(atom, None, label, findings)
            if spec is not None:
                plan[index] = spec
                residue.heavy.append(index)
            continue

        spec = _atom_spec(atom, row, f'{label} {key_name}', findings)
        if spec is None:
            continue
        plan[index] = spec
        residue.by_name[key_name] = index
        residue.heavy.append(index)

    residue.missing = sorted(n for n in residue.template.atoms if n not in residue.by_name)


def _atom_spec(atom, row, where: str, findings: _Findings) -> tuple | None:
    """`(element, charge, isotope)` for one atom, or `None` when it cannot be stored at all.

    The file's element wins -- it is a statement about this atom, the table's about a component -- and
    the table is consulted only where the file stated none.  The charge goes the other way: the record
    cannot tell a stated zero from a blank column, so a zero yields to the table and a non-zero wins.
    """
    element, charge = (row if row is not None else (None, 0))
    if atom.element is None:
        if element is None:
            findings.hit('atom', 'atom(s) state no element and match no template atom, so they cannot '
                                 'be stored at all', where)
            return None
        findings.hit('atom', f"atom(s) state no element; the template's element {element} is used",
                     where)
    else:
        if element is not None and atom.element != element:
            findings.hit('atom', f'atom(s) state element {atom.element} where the template has '
                                 f"{element}; the file's element is used", where)
        element = atom.element

    if atom.charge:
        if atom.charge != charge:
            findings.hit('atom', f'atom(s) state charge {atom.charge} where the template has {charge}; '
                                 f"the file's charge is used", where)
        charge = atom.charge
    return element, charge, atom.isotope


def _plan_untemplated(record: PDBRecord, indexes: list[int], plan: dict[int, tuple],
                      loose: set[int], label: str, findings: _Findings) -> None:
    """Decide what to store for a residue the table has no row for -- a ligand.

    Every atom the file gave is kept, hydrogens included, and only the bonds the file stated are
    applied.  `chython.chemistry.saturate()` is the caller's separate pass and is not called here.
    """
    for index in indexes:
        spec = _atom_spec(record.atoms[index], None, label, findings)
        if spec is not None:
            plan[index] = spec
            loose.add(index)


def _stated_bonds(record: PDBRecord, pairs: dict[tuple[int, int], int],
                  from_template: set[tuple[int, int]], plan: dict[int, tuple], kept: set[int],
                  dropped_h: set[int], findings: _Findings) -> None:
    """Fold the bonds the file stated into the template's, deduplicating by the unordered pair.

    `PDBBond.stated_order` decides three cases: the file stated an order, so the file wins and a
    disagreement is logged; the file named the bond with no order (`CONECT`, `SSBOND`, `LINK`, a
    `_struct_conn` row with no `value_order`, all of which carry order 1 by convention), so the template
    wins silently; or only the template names the bond, and it is applied.  Deduplication by
    `PDBBond.key` is not hypothetical -- `_chem_comp_bond` restates the template's own bonds.
    """
    for bond in record.bonds:
        if bond.a == bond.b:
            continue
        if bond.a in dropped_h or bond.b in dropped_h:
            continue          # our own drop, already reported as a count; not the file's doing
        if bond.a not in kept or bond.b not in kept:
            continue          # our own conformer selection, likewise
        if bond.a not in plan or bond.b not in plan:
            findings.hit('bond', 'stated bond(s) name an atom that could not be stored, so the bond '
                                 'is not built', f'atoms {bond.a}-{bond.b}')
            continue

        key = bond.key
        if key in from_template:
            if bond.stated_order and bond.order != pairs[key]:
                findings.hit('bond', f'stated bond(s) state order {bond.order} where the template has '
                                     f"{pairs[key]}; the file's order is used", f'atoms {key[0]}-{key[1]}')
                pairs[key] = bond.order
        elif key not in pairs:
            pairs[key] = bond.order


def _chain_links(residues: dict[tuple, _Residue], pairs: dict[tuple[int, int], int],
                 findings: _Findings) -> None:
    """Bond `link_out` of each residue to `link_in` of the next one along its chain.

    A residue takes part iff its row names a link atom, which keeps the polymer-kind vocabulary in
    `_residues.py` alone.  The chain is sorted and never taken in file order.  Two residues are joined
    when they are adjacent in that sort, their sequence numbers differ by 1 (or 0, an insertion code),
    their templates are the same kind, both link atoms are present, and no stated bond already joins
    them.  A gap gets a line naming both residues.
    """
    chains: dict[tuple, list[tuple]] = {}
    for key, residue in residues.items():
        if residue.template.link_in is None and residue.template.link_out is None:
            continue
        if key[3] is None:
            findings.hit('residue', 'polymer residue(s) state no sequence number, so they cannot be '
                                    'ordered within their chain and no chain link is built to them',
                         _label(key))
            continue
        chains.setdefault((key[0], key[1]), []).append(key)

    for group in chains.values():
        group.sort(key=lambda k: (k[3], k[4] or ''))
        for first, second in zip(group, group[1:]):
            one, other = residues[first], residues[second]
            where = f'{_label(first)}-{_label(second)}'
            if one.template.kind != other.template.kind:
                findings.hit('residue', f'adjacent residue pair(s) in one chain are a '
                                        f'{one.template.kind} and a {other.template.kind}; the table '
                                        f'states no link between two kinds, so none is built', where)
                continue
            if second[3] - first[3] not in (0, 1):
                findings.hit('residue', 'residue pair(s) are adjacent in their chain but their '
                                        'sequence numbers skip, so no chain link is built between '
                                        'them', where)
                continue
            a = one.by_name.get(one.template.link_out)
            b = other.by_name.get(other.template.link_in)
            if a is None or b is None:
                findings.hit('residue', f'chain link(s) are not built because a link atom the '
                                        f'template names ({one.template.link_out} or '
                                        f'{other.template.link_in}) is absent', where)
                continue
            key = (a, b) if a <= b else (b, a)
            pairs.setdefault(key, 1)     # a LINK or a struct_conn row may state it already
            one.links_built.add(one.template.link_out)
            other.links_built.add(other.template.link_in)


def _shortfall(residues: dict[tuple, '_Residue'], findings: _Findings) -> None:
    """Report template atoms the file never stated, suppressing the ones a built chain link accounts for.

    A mid-chain polymer residue is legitimately short the atom the chain link displaced -- every amino
    acid row names OXT, every nucleotide row OP3.  The rule: for each link atom where a link was built,
    take the `missing` atoms whose *only* template bond attaches to it, and suppress it if there is
    exactly one.  Two or more suppresses none, which is what keeps a genuinely broken residue audible
    (in ALA both ``O`` and ``OXT`` hang off ``C``).  No atom name appears in this function.
    """
    for residue in residues.values():
        if not residue.missing:
            continue
        # Template neighbours: atom_name -> the template names bonded to it.
        nbrs: dict[str, set[str]] = {}
        for a, b, _ in residue.template.bonds:
            nbrs.setdefault(a, set()).add(b)
            nbrs.setdefault(b, set()).add(a)

        suppressed: set[str] = set()
        for link_atom in residue.links_built:
            candidates = [n for n in residue.missing
                          if nbrs.get(n) == {link_atom}]
            if len(candidates) == 1:
                suppressed.add(candidates[0])

        reportable = [n for n in residue.missing if n not in suppressed]
        if reportable:
            name = residue.key[2]
            findings.hit('residue', f'residue(s) of {name} lack the template atom(s) '
                                    f'{", ".join(reportable)}; no bond to an absent atom is applied '
                                    f'and no atom is invented', _label(residue.key))


def build_molecule(record: 'PDBRecord | Sequence[PDBRecord]', *, alt_loc: str | None = None,
                   log: list[str] | None = None) -> MoleculeContainer:
    """One :class:`PDBRecord` -- one model -- as a :class:`~chython.core.MoleculeContainer`.

    A residue the table knows gets the bonds and orders the table states; one it does not know keeps
    every atom the file gave and only the bonds the file stated.  The container is normally disconnected
    (protein, ligands, waters) and a caller who wants the pieces splits it.  `alt_loc` selects an
    alternate conformer: `None` keeps every conformer-free atom plus the conformer with the highest
    summed occupancy, a string keeps that id and logs every residue lacking it.  Atoms of different
    conformers are never bonded.  Three things it does not do: no projection or rotation (x and y go to
    `SEG_XY` as they stand, all three to `SEG_CONFORMERS`; `clean2d()` is the caller's layout pass), no
    distance check on the chain link, and no `saturate()`.

    THE MOLECULE'S OWN `log` IS THE DESTINATION, unconditionally, and it holds the record's parse log
    ahead of this pass's own lines: the reader's findings -- an element column that was not a symbol, a
    CONECT naming an absent serial -- are about these atoms, and `pdb()`'s caller-supplied list is
    someone else's copy.  `log=` here receives this pass's lines only, findings counted by kind
    (`unsupported: ` means our table is the limitation, any other prefix means a broken file).

    A SEQUENCE OF RECORDS COLLAPSES TO ONE MOLECULE WITH A CONFORMER EACH.  Passing a list is the
    opt-in and nothing collapses without it: `pdb()` and `mmcif()` yield one record per model, which is
    what the file said.  `records[0]` builds the molecule and its layout; each further record is matched
    atom for atom against it on `(chain, residue_seq, ins_code, residue_name, atom_name, alt_loc)` and
    becomes one conformer, carrying that record's `model` as the conformer's `ext_index`.  A record
    whose atom set does not match is logged and skipped, the unit of all-or-nothing being one model.
    """
    log = [] if log is None else log
    records = [record] if isinstance(record, PDBRecord) else list(record)
    if not records:
        raise ValueError('build_molecule needs at least one record; an empty sequence states no '
                         'molecule to build')
    record = records[0]
    # This pass's own lines.  Kept apart from `log` so the absorb below is this record's and no other's.
    own: list = []
    findings = _Findings()

    # Imported here, not at module scope, to keep `import chython.formats` cheap.
    from ...chemistry import residue_template

    groups = _residue_groups(record)
    kept = _select_conformers(record, groups, alt_loc, findings)

    plan: dict[int, tuple] = {}                       # record index -> (element, charge, isotope)
    pairs: dict[tuple[int, int], int] = {}            # unordered record-index pair -> bond order
    residues: dict[tuple, _Residue] = {}
    dropped_h: set[int] = set()
    loose: set[int] = set()                           # atoms of residues with no template row
    untemplated: dict[str, int] = {}

    for key, indexes in groups.items():
        selected = [i for i in indexes if i in kept]
        if not selected:
            continue
        name = key[2]
        template = residue_template(name) if name else None
        if template is None:
            untemplated[name or '?'] = untemplated.get(name or '?', 0) + 1
            _plan_untemplated(record, selected, plan, loose, _label(key), findings)
            continue
        residue = _Residue(key, template)
        residues[key] = residue
        _plan_templated(record, residue, selected, plan, dropped_h, findings)
        for one, other, order in template.bonds:
            a, b = residue.by_name.get(one), residue.by_name.get(other)
            if a is None or b is None:
                continue                              # both endpoints or no bond, and no invented atom
            pairs[(a, b) if a <= b else (b, a)] = order

    from_template = set(pairs)
    _stated_bonds(record, pairs, from_template, plan, kept, dropped_h, findings)
    _chain_links(residues, pairs, findings)
    _shortfall(residues, findings)

    mol = MoleculeContainer()
    sids: dict[int, int] = {}
    with mol.edit():
        # Record order, so the container's atom order is the file's.
        for index in sorted(plan):
            element, charge, isotope = plan[index]
            for drop in ((), ('isotope',), ('isotope', 'charge')):
                try:
                    sids[index] = mol.add_atom(element, charge=0 if 'charge' in drop else charge,
                                               isotope=0 if 'isotope' in drop else isotope)
                except ValueError as e:
                    reason = e
                    continue
                if drop:
                    findings.hit('atom', f'atom(s) cannot be stored as stated ({reason}); '
                                         f'{" and ".join(drop)} dropped',
                                 _label(record.atoms[index].residue_key))
                break
            else:
                findings.hit('atom', f'atom(s) cannot be stored even without their isotope and charge '
                                     f'({reason}), so they are not stored',
                             _label(record.atoms[index].residue_key))
        for (a, b), order in pairs.items():
            if a in sids and b in sids:
                mol.add_bond(sids[a], sids[b], order)

    _coordinates(records, mol, sids, findings, alt_loc)
    _hydrogens(record, mol, sids, residues, findings, own)

    findings.report(own)
    if untemplated:
        named = ', '.join(f'{name} ({count})' for name, count in sorted(untemplated.items()))
        own.append(LogRecord('pdb:untemplated-residues', (),
                             f'unsupported: {sum(untemplated.values())} residue(s) of '
                             f'{len(untemplated)} component id(s) have no row in the residue table, '
                             f'so no template bond is applied to them: {named}',
                             LOST))
    unbonded = sum(1 for index in loose if index in sids and not mol.degree_of(sids[index]))
    if unbonded:
        own.append(LogRecord('pdb:unbonded-untemplated', (),
                             f'residue: {unbonded} atom(s) of residue(s) with no template row hold '
                             f'no bond at all; nothing is invented for them, and '
                             f'chython.chemistry.saturate() is the separate pass a caller runs on a '
                             f'ligand whose file gave connectivity and no orders'))
    if record.title or record.entry_id:
        mol.set_title(record.title or record.entry_id)
    # The record's parse log first, then this pass's: the molecule reads in the order it was built.
    mol.log.absorb('read', record.log)
    mol.log.absorb('read', own)
    log.extend(own)
    return mol


def _coordinates(records: list, mol: MoleculeContainer, sids: dict[int, int],
                 findings: _Findings, alt_loc: str | None) -> None:
    """Write the coordinates in a second edit scope, the way every builder in this tree does.

    Both segments: `SEG_XY` is the depiction `clean2d()` may replace, `SEG_CONFORMERS` the stated
    geometry it may not.  An atom with x and y but no z is placed in the plane only, no z invented.  A
    missing x or y is reported only when some other atom has one, a record with no coordinates at all
    having already been reported by the reader that produced it.

    Model 0 is `records[0]`, added explicitly rather than left to the first `set_xyz` so that it carries
    that record's `MODEL` number; a file with no `MODEL` card states None and stores the sentinel.
    """
    record = records[0]
    placed, unplaced = [], []
    for index in sorted(sids):
        atom = record.atoms[index]
        (placed if atom.x is not None and atom.y is not None else unplaced).append(index)
    if not placed:
        return
    # Every z zero means no geometry, so no conformer segment: the same test `_ctab` and `mol2` apply.
    solid = any(record.atoms[index].z for index in placed
                if record.atoms[index].z is not None)
    with mol.edit():
        model = mol.add_conformer(ext_index=_ext_index(record, findings)) if solid else 0
        for index in placed:
            atom = record.atoms[index]
            try:
                mol.set_xy(sids[index], atom.x, atom.y)
                if solid and atom.z is not None:
                    mol.set_xyz(sids[index], atom.x, atom.y, atom.z, model=model)
            except ValueError as e:
                findings.hit('coordinates', f'atom(s) hold a coordinate the container refuses ({e}), '
                                            f'so none is stored for them',
                             _label(atom.residue_key))
    for index in unplaced:
        findings.hit('coordinates', 'atom(s) state no x or y coordinate where other atoms of this '
                                    'record do, so no position is stored for them',
                     _label(record.atoms[index].residue_key))
    if solid:
        # A conformer segment is one dense column per model with no per-atom validity flag, so an atom
        # with no z reads back at the origin, indistinguishable from one the file put there.  Reported
        # here because the segment cannot report it.
        for index in placed:
            if record.atoms[index].z is None:
                findings.hit('coordinates', 'atom(s) state x and y but no z where other atoms of this '
                                            'record state one, so they sit at the origin in the '
                                            'stored geometry',
                             _label(record.atoms[index].residue_key))
    if solid and len(records) > 1:
        _extra_models(records, mol, sids, placed, alt_loc, findings)


def _extra_models(records: list, mol: MoleculeContainer, sids: dict[int, int], placed: list,
                  alt_loc: str | None, findings: _Findings) -> None:
    """Every record after the first as one further conformer, atoms matched by annotation.

    ALL-OR-NOTHING PER MODEL: a record whose kept atom set does not match `records[0]`'s is logged and
    skipped, and the rest still land.  The edit scope is what enforces it for a coordinate the container
    refuses -- leaving the `with` by exception discards the journal, so a half-filled model cannot be
    sealed.
    """
    # The same filter on both sides -- a placed atom with no z sits at the origin in model 0 and is not
    # part of the set being matched, so one such atom does not disqualify every further model.
    keys = {_match_key(records[0].atoms[index]): index
            for index in placed if records[0].atoms[index].z is not None}
    for record in records[1:]:
        groups = _residue_groups(record)
        kept = _select_conformers(record, groups, alt_loc, findings)
        by_key = {}
        for index in sorted(kept):
            atom = record.atoms[index]
            if atom.x is None or atom.y is None or atom.z is None:
                continue
            by_key[_match_key(atom)] = index
        if set(by_key) != set(keys):
            findings.hit('conformer', 'model(s) state a different atom set from the first model, so no '
                                      'conformer is stored for them', _model_label(record))
            continue
        try:
            with mol.edit():
                model = mol.add_conformer(ext_index=_ext_index(record, findings))
                for key, index in by_key.items():
                    atom = record.atoms[index]
                    mol.set_xyz(sids[keys[key]], atom.x, atom.y, atom.z, model=model)
        except ValueError as e:
            findings.hit('conformer', f'model(s) hold a coordinate the container refuses ({e}), so no '
                                      f'conformer is stored for them', _model_label(record))


def _hydrogens(record: PDBRecord, mol: MoleculeContainer, sids: dict[int, int],
               residues: dict[tuple, _Residue], findings: _Findings, log: list[str]) -> None:
    """Derive the implicit hydrogen counts, then check the drop against the derivation.

    The one shared derivation, run outside any edit scope because it needs a sealed arena, filling only
    what nothing has claimed.  A templated residue's explicit hydrogens were dropped, so where the stated
    count differs from the derived one the difference is reported, aggregated by residue name and the two
    counts.  A residue the file gave no explicit hydrogens for is not checked: a heavy-atom-only file --
    almost every archive entry -- states nothing about protonation.
    """
    unsettled = mol.derive_hydrogens()
    if unsettled:
        log.append(LogRecord('pdb:unsettled-hydrogens', tuple(unsettled),
                             f'atom: {len(unsettled)} atom(s) hold no derivable implicit hydrogen '
                             f'count; kekule() settles the aromatic pnictogen and check_valence() '
                             f'names the rest',
                             LOST))
    for key, residue in residues.items():
        if not residue.stated_h:
            continue
        # One unsettled count would understate the sum below and fire the finding falsely.
        if any(mol.implicit_h_of(sids[i]) is None for i in residue.heavy if i in sids):
            continue
        derived = sum(mol.implicit_h_of(sids[i]) for i in residue.heavy if i in sids)
        if derived != residue.stated_h:
            findings.hit('atom', f'residue(s) of {key[2]} state {residue.stated_h} explicit '
                                 f'hydrogen(s) where {derived} are derived; the explicit hydrogens '
                                 f'are dropped and the derived count is used', _label(key))
