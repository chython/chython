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
"""The neutral intermediate both PDB-family readers produce: atoms, stated bonds, annotation, log.

A field the file did not state is ``None``, never a zero and never an empty string -- an occupancy of
0.0 is a statement about a disordered atom, an omitted one is not.  Two exceptions: ``charge`` is 0 and
``hetatm`` is False, so an mmCIF with no ``group_PDB`` column (legal) reads as all-ATOM and cannot be
told from one that stated ``ATOM`` on every row; ``entity_type`` is the field that stays ``None``.
Of the four author identifiers only ``auth_chain``/``auth_seq`` are kept -- the numbering is what
differs from the label spelling in practice -- and ``auth_comp_id``/``auth_atom_id`` are logged unread.
"""

from ...core._core import element_symbols


__all__ = ['PDBAtom', 'PDBBond', 'PDBRecord', 'normalize_element', 'range_messages']


_SYMBOLS = element_symbols()                      # ('R', 'H', 'He', ..., 'Og')
_VALID = frozenset(_SYMBOLS[1:])
_UPPER = {s.upper(): s for s in _VALID}           # 'FE' -> 'Fe', 'CL' -> 'Cl'


def normalize_element(raw: str) -> tuple[str | None, int, str | None]:
    """``(symbol, isotope, log_message_or_None)`` for an element field from either PDB dialect.

    Both dialects write the symbol in upper case -- legacy PDB in a right-justified two-column field,
    mmCIF in ``_atom_site.type_symbol`` -- so the case fold is not optional.  An unrecognised field
    yields ``None`` and a message, never a guess from the atom name: ``NA`` in the atom-name column of a
    haem is a pyrrole nitrogen, not sodium.
    """
    token = raw.strip()
    if not token:
        return None, 0, None                      # the caller knows whether an empty field is news
    if token in _VALID:
        return token, 0, None
    upper = token.upper()
    symbol = _UPPER.get(upper)
    if symbol is not None:
        return symbol, 0, None                    # the expected spelling in both formats, not damage
    if upper == 'D':
        return 'H', 2, "atom: element 'D' read as hydrogen, isotope 2"
    if upper == 'T':
        return 'H', 3, "atom: element 'T' read as hydrogen, isotope 3"
    # A core-CIF style `Fe2+` or a writer's `C 1`: keep the leading alphabetic run and say so.
    head = ''
    for ch in token:
        if ch.isalpha():
            head += ch
        else:
            break
    if head:
        symbol = _UPPER.get(head.upper())
        if symbol is not None:
            return symbol, 0, (f'atom: element field {token!r} carries more than a symbol; read as '
                               f'{symbol!r}')
    return None, 0, f'atom: element field {token!r} is not an element symbol; element not stated'


def range_messages(occupancy: float | None, b_factor: float | None, where: str) \
        -> list[tuple[str, str]]:
    """``(kind, line)`` pairs for an occupancy or a B factor outside the range the quantity can have.

    Both values are stored as stated; this only reports them.  *kind* is the aggregation key both readers
    count on.  An occupancy is a fraction (0..1) and a negative B factor is wrong in both dialects; a B
    factor over 999.99 is wrong only in legacy PDB, where it does not fit the six-column field, so the
    legacy reader reports that one and this shared check does not.
    """
    messages = []
    if occupancy is not None and not 0.0 <= occupancy <= 1.0:
        messages.append(('an occupancy outside 0.0-1.0',
                         f'atom: occupancy {occupancy} on {where} is outside 0.0-1.0; stored as '
                         f'stated'))
    if b_factor is not None and b_factor < 0.0:
        messages.append(('a negative B factor',
                         f'atom: B factor {b_factor} on {where} is negative; stored as stated'))
    return messages


class PDBAtom:
    """One atom as the file stated it, with its annotation.

    ``element`` is chython's spelling (``'Fe'``, not ``'FE'``) or ``None``, never derived from the atom
    name; ``x``, ``y``, ``z`` are Angstroms, each ``None`` if missing or unreadable.  ``chain`` is
    ``label_asym_id`` and is a string, not a character -- a one-character assumption breaks on large
    assemblies.  ``auth_chain``/``auth_seq`` are the depositor's numbering, which routinely differs from
    the label numbering that ``residue_name`` and ``atom_name`` always carry.  ``alt_loc`` is the
    alternate-conformer id, kept because no reader may bond across conformers.  ``entity_type`` is
    mmCIF's vocabulary (``'polymer'``, ``'non-polymer'``, ``'water'``, ``'branched'``), of which
    ``hetatm`` is a separate statement.
    """
    __slots__ = ('element', 'isotope', 'charge', 'x', 'y', 'z', 'serial', 'atom_name',
                 'residue_name', 'chain', 'auth_chain', 'auth_seq', 'residue_seq', 'ins_code',
                 'alt_loc', 'occupancy', 'b_factor', 'hetatm', 'entity_type', 'model')

    def __init__(self, element: str | None = None, x: float | None = None,
                 y: float | None = None, z: float | None = None, *, isotope: int = 0,
                 charge: int = 0,
                 serial: int | None = None, atom_name: str | None = None,
                 residue_name: str | None = None, chain: str | None = None,
                 auth_chain: str | None = None, auth_seq: int | None = None,
                 residue_seq: int | None = None, ins_code: str | None = None,
                 alt_loc: str | None = None, occupancy: float | None = None,
                 b_factor: float | None = None, hetatm: bool = False,
                 entity_type: str | None = None, model: int | None = None):
        self.element = element
        self.isotope = isotope
        self.charge = charge
        self.x = x
        self.y = y
        self.z = z
        self.serial = serial
        self.atom_name = atom_name
        self.residue_name = residue_name
        self.chain = chain
        self.auth_chain = auth_chain
        self.auth_seq = auth_seq
        self.residue_seq = residue_seq
        self.ins_code = ins_code
        self.alt_loc = alt_loc
        self.occupancy = occupancy
        self.b_factor = b_factor
        self.hetatm = hetatm
        self.entity_type = entity_type
        self.model = model

    @property
    def is_water(self) -> bool:
        """True only when the file stated this atom belongs to a water entity."""
        return self.entity_type == 'water'

    @property
    def is_polymer(self) -> bool:
        """True only when the file stated this atom belongs to a polymer entity."""
        return self.entity_type == 'polymer'

    @property
    def is_ligand(self) -> bool:
        """True only when the file stated a non-polymer, non-water entity."""
        return self.entity_type in ('non-polymer', 'branched')

    @property
    def residue_key(self) -> tuple:
        """The tuple two atoms of one residue share, and two residues never do.

        ``alt_loc`` is deliberately not part of it: alternate conformers describe one residue.  Not
        bonding across them is a rule about bonds and lives in the readers.
        """
        return self.model, self.chain, self.residue_name, self.residue_seq, self.ins_code

    def __repr__(self):
        return (f'PDBAtom({self.element!r}, {self.x!r}, {self.y!r}, {self.z!r}, '
                f'atom_name={self.atom_name!r}, residue_name={self.residue_name!r}, '
                f'chain={self.chain!r}, residue_seq={self.residue_seq!r}, '
                f'alt_loc={self.alt_loc!r})')


class PDBBond:
    """One bond the file stated, never one a reader inferred.

    ``a`` and ``b`` are indices into :attr:`PDBRecord.atoms`.  ``order`` is chython's bond order (1, 2,
    3, 4 aromatic, 8 dative).  ``stated_order`` is False when the file named the bond but no order --
    ``CONECT``, ``SSBOND``, ``LINK`` -- where order 1 is stored and logged and a repeated ``CONECT`` pair
    is never read as a multiplicity, writers disagreeing about what that means.  ``source`` names the
    table or record (``'chem_comp_bond'``, ``'struct_conn'``, ``'conect'``, ``'ssbond'``, ``'link'``),
    ``conn_type`` is ``_struct_conn.conn_type_id`` verbatim (``'disulf'``, ``'covale'``, ``'metalc'``).
    """
    __slots__ = ('a', 'b', 'order', 'stated_order', 'source', 'conn_type')

    def __init__(self, a: int, b: int, order: int = 1, *, stated_order: bool = True,
                 source: str = '', conn_type: str | None = None):
        self.a = a
        self.b = b
        self.order = order
        self.stated_order = stated_order
        self.source = source
        self.conn_type = conn_type

    @property
    def key(self) -> tuple:
        """The unordered atom pair, low index first -- the identity a duplicate is detected by."""
        return (self.a, self.b) if self.a <= self.b else (self.b, self.a)

    def __repr__(self):
        return (f'PDBBond({self.a}, {self.b}, {self.order}, stated_order={self.stated_order}, '
                f'source={self.source!r}, conn_type={self.conn_type!r})')


class PDBRecord:
    """One model out of one file: its atoms, the bonds the file stated, and the parse log.

    A multi-``MODEL`` legacy PDB and a multi-``pdbx_PDB_model_num`` mmCIF both yield one record per
    model -- an NMR ensemble is alternative structures, not one structure with duplicated atoms.
    ``entry_id`` is the ``data_`` block name or the legacy ``HEADER`` id, ``title`` the entry title.
    ``log`` is this record's messages, also appended to the reader's ``log=`` list when one was passed.
    It stays on the record because a record is not a container: :func:`build_molecule` folds it onto the
    ``log`` of the molecule it builds, which is where a caller holding only the molecule reads it.
    """
    __slots__ = ('atoms', 'bonds', 'log', 'entry_id', 'title', 'model')

    def __init__(self, *, entry_id: str | None = None, title: str | None = None,
                 model: int | None = None):
        self.atoms: list[PDBAtom] = []
        self.bonds: list[PDBBond] = []
        self.log: list[str] = []
        self.entry_id = entry_id
        self.title = title
        self.model = model

    def __len__(self):
        return len(self.atoms)

    def unbonded_count(self) -> int:
        """How many atoms no stated bond touches.

        A method and not a cached property: readers call it while still appending bonds.
        """
        bonded = set()
        for bond in self.bonds:
            bonded.add(bond.a)
            bonded.add(bond.b)
        return len(self.atoms) - len(bonded)

    def __repr__(self):
        return (f'PDBRecord(entry_id={self.entry_id!r}, model={self.model!r}, '
                f'atoms={len(self.atoms)}, bonds={len(self.bonds)})')
