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
"""The monomer library: Pistoia JSON records, the process-wide set, and its compiled index.

A record's ``smiles`` marks each **leaving atom** with its R number as a map number -- ``N[H:1]``,
``C([OH:2])=O`` -- so the record states both the free monomer and where it attaches.  The set in use is
the shipped public library (Pistoia HELMCoreLibrary, MIT) plus whatever :func:`set_monomer_library`
attached, later files overriding earlier ones by ``(polymer type, symbol)``.  Nothing loads at import.

| Damage in a record | What loading does |
|---|---|
| SMILES that does not parse, or no mapped atom | record left out, ``polymers:monomer-unreadable`` |
| a symbol already defined in the same file | first kept, ``polymers:monomer-duplicate`` |
"""

from collections import Counter
from gzip import open as gzip_open
from itertools import combinations
from re import compile as re_compile
from json import load as json_load
from pathlib import Path
from typing import NamedTuple

from ._shape import shape_of_core
from ..core import LogRecord, LOST, smiles


__all__ = ['Monomer', 'MonomerLibrary', 'load_library', 'set_monomer_library', 'get_monomer_library',
           'DEFAULT_LIBRARY', 'piece_keys', 'signature', 'markers_for_closures']


#: The shipped public set.
DEFAULT_LIBRARY = Path(__file__).parent / 'data' / 'HELMCoreLibrary.json.gz'

_LIBRARY_CACHE = {}
_ATTACHED = []


class Monomer(NamedTuple):
    symbol: str
    polymer_type: str
    monomer_type: str
    smiles: str
    labels: tuple          # R numbers, ascending
    natural_analog: str


class Fragment(NamedTuple):
    """A monomer's SMILES with each used R written as the token ``{R<n>}`` and its leaving group gone.

    An unused R keeps its leaving group -- the cap.  `closures` are the ring numbers the text uses for its
    own rings, which an assembly must not hand out as attachment ids; `atom_count` is the atoms it writes.
    """
    text: str
    closures: tuple
    atom_count: int


class MonomerLibrary:
    """The monomers of one library set, by ``(polymer type, symbol)``, and the load log."""
    __slots__ = ('_by_symbol', 'monomers', 'log', '_parts', '_fragments', '_index')

    def __init__(self, monomers, log):
        self._by_symbol = {(m.polymer_type, m.symbol): m for m in monomers}
        self.monomers = tuple(self._by_symbol.values())
        self.log = log
        self._parts = {}
        self._fragments = {}
        self._index = None

    def get(self, polymer_type, symbol):
        """The monomer, or ``None``."""
        return self._by_symbol.get((polymer_type, symbol))

    def by_core(self, key):
        """``((monomer, used R numbers), ...)`` whose core has this key, stereo in.

        A core is the monomer with each USED R's leaving group replaced by an unlabelled marker and every
        unused R keeping its cap -- so a terminal residue, cap and all, has a key too.
        """
        return self._built().core.get(key, ())

    def by_shape(self, key):
        """``((monomer, used R numbers), ...)`` whose core has this :mod:`._shape` key -- a screen, not a match."""
        return self._built().shapes.get(key, ())

    def by_constitution(self, key):
        """As :meth:`by_core`, with stereo removed from the key."""
        return self._built().constitution.get(key, ())

    @property
    def core_formulas(self):
        """Every core's element counts, markers as ``R``: what a unit must have to be worth a lookup."""
        return self._built().formulas

    @property
    def cut_scopes(self):
        """``((polymer type, pairs), ...)``: the bond-end pairs a decomposition cuts, one set per type.

        A pair is two attachment signatures of monomers with two or more R that can join: an R2 and an R1
        of the type (the in-chain link), an R3-or-higher of the type with any of its R (a cross-link or a
        branch), or a CHEM monomer's R with any R of the type.  A one-R monomer is a terminal cap, whose
        attachment is an ordinary CH3, CH2 or NH2: it cuts by :meth:`caps` instead.  One end looking like
        some attachment is not enough -- every CH2-CH2 bond of a lysine would be cut -- and a type's small
        linkers must not cut another type's side chains, which is why the sets are per type.  The last scope
        is every type at once.
        """
        return self._built().pairs

    def caps_by_shape(self, key):
        """``((monomer, used), ...)`` for a terminal cap whose core has this :mod:`._shape` key."""
        return self._built().cap_shapes.get(key, ())

    def caps(self, key):
        """``((monomer, used), ...)`` for a one-R monomer -- a terminal cap -- whose core has this key."""
        return self._built().caps.get(key, ())

    @property
    def cap_sites(self):
        """Every R1 and R2 link-end signature: the near end a cap cut must have, checked before anything."""
        built = self._built()
        return frozenset(s for (_, kind), sigs in built.ends.items() if kind in ('R1', 'R2') for s in sigs)

    def link_ends(self, polymer_type, kind):
        """The signatures of `kind` ('R1', 'R2' or 'R3+') attachment atoms of the type's linking monomers."""
        return self._built().ends.get((polymer_type, kind), frozenset())

    @property
    def cap_formulas(self):
        """Element counts of every cap core, markers as ``R``."""
        return self._built().cap_formulas

    @property
    def cut_signatures(self):
        """The :func:`signature` of every attachment atom in the library: what a molecule bond is cut by."""
        return self._built().signatures

    def labelled_core(self, m, used):
        """``(core, {marker atom: R number})``: the core of `m` with `used` open, its markers named."""
        return self._built().labelled[(m.polymer_type, m.symbol, used)]

    def _built(self):
        if self._index is None:
            self._index = _build_index(self)
        return self._index

    def parts(self, m):
        """``(molecule with map numbers cleared, {R number: (attachment atom, leaving atom)})``."""
        key = (m.polymer_type, m.symbol)
        parts = self._parts.get(key)
        if parts is None:
            parts = self._parts[key] = _split_monomer(m.smiles)
        return parts

    def fragment(self, m, used):
        """The :class:`Fragment` of monomer `m` with the R numbers in `used` open; cached."""
        key = (m.polymer_type, m.symbol, used)
        fragment = self._fragments.get(key)
        if fragment is None:
            mol, sites = self.parts(m)
            labels = sorted(used)
            detached = mol.detached_smiles({10 + i: sites[label] for i, label in enumerate(labels)})
            text = detached.text
            for i, label in enumerate(labels):
                text = text.replace(f'%{10 + i}', f'{{R{label}}}')
            fragment = self._fragments[key] = Fragment(text, detached.closure_ids, detached.atom_count)
        return fragment


def _split_monomer(text):
    """The monomer with map numbers cleared, and each R's attachment and leaving atoms.

    The leaving atom is the one carrying the R number.  Its attachment is the neighbour across the bond
    that cuts off a side holding no other leaving atom, the smallest such side when there are several --
    a one-atom cap like ``[OH:2]`` has one neighbour, an azide cap's first nitrogen has two.
    """
    mol = smiles(text)
    leaving = {mol.atom(n).map_number: n for n in mol.atom_numbers if mol.atom(n).map_number}
    sites = {}
    for label, atom in leaving.items():
        best = None
        for neighbour in mol.neighbors_of(atom):
            side = _side(mol, atom, neighbour)
            if neighbour in side or any(a in side for a in leaving.values() if a != atom):
                continue
            if best is None or len(side) < best[0]:
                best = (len(side), neighbour)
        if best is None:
            raise ValueError(f'R{label} of {text!r} has no attachment atom')
        sites[label] = (best[1], atom)
    with mol.edit():
        for atom in leaving.values():
            mol.set_map_number(atom, 0)
    return mol, sites


def _side(mol, start, across):
    """The atoms reachable from `start` without crossing the bond to `across`."""
    seen, stack = {start}, [start]
    while stack:
        atom = stack.pop()
        for n in mol.neighbors_of(atom):
            if n not in seen and not (atom == start and n == across):
                seen.add(n)
                stack.append(n)
    return seen


def load_library(path):
    """``(monomers, log)`` of one Pistoia JSON file, ``.json`` or ``.json.gz``."""
    path = Path(path)
    opener = gzip_open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as f:
        records = json_load(f)
    monomers, log, seen = [], [], set()
    for record in records:
        symbol = str(record.get('symbol') or '')
        polymer_type = str(record.get('polymerType') or '').upper()
        text = record.get('smiles') or ''
        labels = _labels(text)
        if not symbol or not labels:
            log.append(LogRecord('polymers:monomer-unreadable', (),
                                 f'{path.name}: monomer {symbol or "?"} ({polymer_type or "?"}): SMILES {text!r:.60} '
                                 f'does not parse or marks no leaving atom; left out', LOST))
            continue
        key = (polymer_type, symbol)
        if key in seen:
            log.append(LogRecord('polymers:monomer-duplicate', (),
                                 f'{path.name}: monomer {symbol} ({polymer_type}) defined again; the first record '
                                 f'is kept', LOST))
            continue
        seen.add(key)
        monomers.append(Monomer(symbol, polymer_type, str(record.get('monomerType') or ''), text, labels,
                                str(record.get('naturalAnalog') or '')))
    return monomers, log


def _labels(text):
    """The R numbers a record's SMILES marks, ascending; empty when it does not parse."""
    try:
        mol = smiles(text)
    except (ValueError, TypeError):
        return ()
    return tuple(sorted(mol.atom(n).map_number for n in mol.atom_numbers if mol.atom(n).map_number))


def set_monomer_library(*paths):
    """Use the shipped library plus `paths`, in order, later overriding earlier.  No paths: the default."""
    _ATTACHED[:] = [Path(p) for p in paths]
    _LIBRARY_CACHE.clear()


def get_monomer_library():
    """The library set in use, built once and cached until :func:`set_monomer_library` changes it."""
    library = _LIBRARY_CACHE.get('set')
    if library is None:
        merged, log = {}, []
        for path in (DEFAULT_LIBRARY, *_ATTACHED):
            monomers, own = load_library(path)
            log.extend(own)
            for m in monomers:
                merged.pop((m.polymer_type, m.symbol), None)     # an override moves to the end: later wins
                merged[(m.polymer_type, m.symbol)] = m
        library = _LIBRARY_CACHE['set'] = MonomerLibrary(list(merged.values()), log)
    return library


_SCOPE_ORDER = {'PEPTIDE': 0, 'RNA': 1, 'DNA': 2, 'CHEM': 3}
_CLOSURE = re_compile(r'([-=#:/\\]?)%([0-9]{2})')


class _Index(NamedTuple):
    core: dict
    constitution: dict
    signatures: frozenset
    labelled: dict
    formulas: frozenset
    pairs: tuple
    caps: dict
    cap_formulas: frozenset
    ends: dict
    shapes: dict
    cap_shapes: dict


def markers_for_closures(text):
    """`text` with every two-digit closure ``%nn`` written as a ``(*)`` branch, and the ids in text order.

    A ring-bond digit and a branch take the same place in an atom's neighbour order, so a configuration
    written around the closure reads the same around the marker.
    """
    order = []

    def branch(match):
        order.append(int(match.group(2)))
        return f'({match.group(1)}*)'

    return _CLOSURE.sub(branch, text), order


def piece_keys(text):
    """``(key, constitution key)`` of a SMILES with ``*`` markers: canonical bytes with and without stereo.

    Aromaticity is normalised first (`thiele`), so a Kekule and an aromatic drawing of one piece share a key.
    """
    mol = smiles(text)
    return _keys(mol)


def _keys(mol):
    mol.implicify_hydrogens()       # a kept cap's hydrogen, written explicitly by the library record
    mol.thiele()
    flat = mol.copy()
    flat.clean_stereo()
    return mol.canonical_bytes, flat.canonical_bytes


def signature(mol, atom, other):
    """What one end of the bond `atom`-`other` looks like with that bond removed.

    ``(element, aromatic, charge, implicit hydrogens, heavy neighbours other than `other`, heteroatom
    neighbours other than `other`)``.  The last field is what tells a triazole carbon, beside its ring
    nitrogen, from a benzene ipso carbon -- without it every benzyl bond reads as a click junction.
    """
    a = mol.atom(atom)
    neighbours = [n for n in mol.neighbors_of(atom) if n != other]
    aromatic = any(mol.order_of(atom, n) == 4 for n in mol.neighbors_of(atom))
    hetero = sum(1 for n in neighbours if mol.atom(n).atomic_symbol not in ('C', 'H', 'R'))
    return a.atomic_symbol, aromatic, a.charge, mol.implicit_h_of(atom), len(neighbours), hetero


def formulas_of(piece):
    """Element counts of a piece, markers counted as ``R``."""
    return tuple(sorted(Counter(piece.atom(n).atomic_symbol for n in piece.atom_numbers).items()))


def _scopes(ends, signatures):
    """One pair set per polymer type, PEPTIDE first, CHEM folded into each; then all types at once."""
    types = sorted({t for t, _ in ends if t != 'CHEM'}, key=lambda t: (_SCOPE_ORDER.get(t, 9), t))
    chem = set()
    for (t, _), sigs in ends.items():
        if t == 'CHEM':
            chem |= sigs
    scopes, everything = [], set()
    for t in types:
        own = set()
        for (u, _), sigs in ends.items():
            if u == t:
                own |= sigs
        pairs = {(a, b) for a in ends.get((t, 'R2'), ()) for b in ends.get((t, 'R1'), ())}
        pairs |= {(a, b) for a in ends.get((t, 'R3+'), ()) for b in own}
        pairs |= {(a, b) for a in chem for b in own}
        pairs |= {(b, a) for a, b in pairs}
        everything |= pairs
        scopes.append((t, frozenset(pairs)))
    if len(scopes) > 1:
        scopes.append(('ALL', frozenset(everything)))
    return scopes


def _build_index(library):
    core, constitution, labelled, signatures, formulas = {}, {}, {}, set(), set()
    caps, cap_formulas = {}, set()
    shapes, cap_shapes = {}, {}
    ends = {}                       # (polymer type, 'R1' | 'R2' | 'R3+') -> signatures
    for m in library.monomers:
        try:
            mol, sites = library.parts(m)
        except ValueError as e:
            library.log.append(LogRecord('polymers:monomer-unindexable', (),
                                         f'monomer {m.symbol} ({m.polymer_type}): {e}; not used to read molecules',
                                         LOST))
            continue
        labels = sorted(sites)
        for size in range(len(labels) + 1):
            for used in combinations(labels, size):
                detached = mol.detached_smiles({10 + i: sites[label] for i, label in enumerate(used)})
                text, order = markers_for_closures(detached.text)
                piece = smiles(text)
                markers = [n for n in piece.atom_numbers if piece.atom(n).is_r]
                names = {atom: used[number - 10] for atom, number in zip(markers, order)}
                key, flat = _keys(piece)
                formulas.add(formulas_of(piece))
                entry = (m, frozenset(used))
                core.setdefault(key, []).append(entry)
                constitution.setdefault(flat, []).append(entry)
                labelled[(m.polymer_type, m.symbol, frozenset(used))] = (piece, names)
                shapes.setdefault(shape_of_core(piece), []).append(entry)
                if size == len(labels) == 1:
                    caps.setdefault(key, []).append(entry)
                    cap_shapes.setdefault(shape_of_core(piece), []).append(entry)
                    cap_formulas.add(formulas_of(piece))
                if size == len(labels) > 1:
                    for marker in markers:
                        (attach,) = piece.neighbors_of(marker)
                        end = signature(piece, attach, marker)
                        signatures.add(end)
                        label = names[marker]
                        kind = 'R1' if label == 1 else 'R2' if label == 2 else 'R3+'
                        ends.setdefault((m.polymer_type, kind), set()).add(end)
    pairs = tuple(_scopes(ends, signatures))
    return _Index({k: tuple(v) for k, v in core.items()}, {k: tuple(v) for k, v in constitution.items()},
                  frozenset(signatures), labelled, frozenset(formulas), pairs,
                  {k: tuple(v) for k, v in caps.items()}, frozenset(cap_formulas),
                  {k: frozenset(v) for k, v in ends.items()}, {k: tuple(v) for k, v in shapes.items()},
                  {k: tuple(v) for k, v in cap_shapes.items()})
