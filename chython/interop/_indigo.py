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
"""
Indigo conversion, both directions.  Reached through `interop.indigo`.

Indigo derives cis/trans configuration from 2D coordinates via `markStereobonds`, so without a layout
cis/trans is a logged loss; its bond-stereo values (7 = Z, 8 = E) likewise cannot be decoded back.
Allene, atropisomer and helical stereo have no Indigo form.  Dative bonds survive in the Indigo
object but not through its SMILES.
"""
from ..core import LogRecord, LOST, SU_CIS_TRANS, SU_TETRA
from ._records import deliver
from ._stereo import kind_name, set_parity_by_probe


def to_indigo(mol, /, *, log=None):
    """
    Export a chython container as an Indigo molecule.

    Atoms are emitted in the molecule's own iteration order and `iterateAtoms()` follows it; that is
    part of the interface, since `depict/layout/molecule.py` maps 2D coordinates back by Indigo atom
    position.  Nothing kekulizes: stored bond orders go out as-is, so call ``thiele()`` first if you
    want aromatic output.  Reportable losses: ``H_UNKNOWN`` (Indigo cannot state it, so it perceives
    a count from valence rules), cis/trans and allene stereo without a 2D layout, and dative bonds
    (order 8), which the Indigo object holds but its SMILES does not.

    :param log: optional list; one human-readable line is appended per reportable loss.
    """
    from indigo import Indigo as _Indigo
    from ..core import MoleculeContainer as _CoreMC

    ig = _Indigo()
    ig_mol = ig.createMolecule()

    is_v3 = isinstance(mol, _CoreMC)

    if is_v3:
        _syms = _CoreMC.__module__  # import side-effect: ensure core is loaded
        from ..core._core import element_symbols as _element_symbols
        _syms = _element_symbols()

        id_to_idx = {}
        has_dative = False
        has_h_unknown = False
        has_r = False

        for a in mol.atoms():
            if a.element == 0:
                has_r = True
            ia = ig_mol.addAtom(a.atomic_symbol)
            if a.charge:
                ia.setCharge(a.charge)
            if a.isotope:
                ia.setIsotope(a.isotope)
            if a.is_radical:
                ia.setRadical(102)  # Indigo: 102 = doublet (mono-radical)
            if a.implicit_h is None:
                # Indigo cannot spell H_UNKNOWN; log and leave the count unset
                has_h_unknown = True
                if log is not None:
                    log.append(LogRecord('indigo:h-unknown-not-representable', (a.n,),
                                         f'atom {a.n} ({a.atomic_symbol}): H_UNKNOWN cannot be stated '
                                         f'in Indigo; hydrogen count not set and will be perceived from valence rules',
                                         LOST))
            else:
                ia.setImplicitHCount(a.implicit_h)
            id_to_idx[a.n] = ia.index()

        for b in mol.bonds():
            if b.order == 8:
                has_dative = True
            ig_mol.getAtom(id_to_idx[b.n]).addBond(ig_mol.getAtom(id_to_idx[b.m]), b.order)

        if has_dative and log is not None:
            log.append(LogRecord(
                'indigo:dative-bonds-lost-in-smiles', (),
                'dative bonds (order 8) are held by the Indigo object but are lost when it is serialised to '
                'SMILES; round-tripping through SMILES will convert them to single bonds', LOST))

        if has_r and log is not None:
            log.append(LogRecord('indigo:r-marker-not-read-back', (),
                                 'an R marker exports as an Indigo pseudo-atom carrying its label; '
                                 '`from_indigo` raises on a pseudo-atom, so this export is one-way',
                                 LOST))

        for u in mol.stereogenic_units():
            if u['kind'] != SU_TETRA or not u['stereogenic'] or u['parity'] == 0:
                continue
            anchor = u['anchor']
            refs = u['refs']
            parity = u['parity']

            # Indigo pyramid: stable ids -> Indigo atom indices; None -> -1 (implicit H)
            pyramid = [id_to_idx[r] if r is not None else -1 for r in refs]

            # Indigo reads the pyramid CW from v1 (ABS type 1), so parity 2 (odd) needs one swap.
            if parity == 2:
                pyramid[0], pyramid[1] = pyramid[1], pyramid[0]

            ig_mol.getAtom(id_to_idx[anchor]).addStereocenter(
                _Indigo.ABS, pyramid[0], pyramid[1], pyramid[2], pyramid[3]
            )

        # Indigo derives cis/trans from 2D geometry, so hand it the layout and let markStereobonds read it
        if mol.has_coordinates:
            for a in mol.atoms():
                xy = a.xy
                if xy is not None:
                    ig_mol.getAtom(id_to_idx[a.n]).setXYZ(xy[0], xy[1], 0.0)
            ig_mol.markStereobonds()
        elif log is not None and any(
            u['kind'] == SU_CIS_TRANS and u['stereogenic'] and u['parity'] != 0
            for u in mol.stereo_units()
        ):
            log.append(LogRecord('indigo:cis-trans-requires-coordinates', (),
                                 'cis/trans stereo cannot be exported without 2D coordinates; '
                                 'call mol.calculate2d() (or another 2D layout method) before exporting, '
                                 'then call markStereobonds() on the Indigo object',
                                 LOST))

        # Every other stereo kind: no Indigo form at all, layout or not -- a layout only lets Indigo
        # read double-bond geometry.  Named from the core's table so a new kind cannot go unreported.
        if log is not None:
            for u in mol.stereo_units():
                kind = u['kind']
                if kind in (SU_TETRA, SU_CIS_TRANS) or not u['stereogenic'] or u['parity'] == 0:
                    continue
                log.append(LogRecord('indigo:stereo-kind-not-representable', (u["anchor"],),
                                     f'atom {u["anchor"]}: {kind_name(kind)} stereo has no Indigo representation; '
                                     f'not exported', LOST))

    else:
        # a chython 2 MoleculeContainer
        mapping = {}
        has_dative = False

        for n, a in mol.atoms():
            ia = ig_mol.addAtom(a.atomic_symbol)
            if a.charge:
                ia.setCharge(a.charge)
            if a.isotope:
                ia.setIsotope(a.isotope)
            if a.is_radical:
                ia.setRadical(102)  # Indigo: 102 = doublet (mono-radical)
            if a.implicit_hydrogens is None:
                if log is not None:
                    log.append(LogRecord('indigo:h-unknown-not-representable', (),
                                         f'atom {n} ({a.atomic_symbol}): H_UNKNOWN cannot be stated in Indigo; '
                                         f'hydrogen count not set and will be perceived from valence rules',
                                         LOST))
            else:
                ia.setImplicitHCount(a.implicit_hydrogens)
            mapping[n] = ia.index()

        for n, m, b in mol.bonds():
            if b.order == 8:
                has_dative = True
            ig_mol.getAtom(mapping[n]).addBond(ig_mol.getAtom(mapping[m]), b.order)

        if has_dative and log is not None:
            log.append(LogRecord(
                'indigo:dative-bonds-lost-in-smiles', (),
                'dative bonds (order 8) are held by the Indigo object but are lost when it is serialised to '
                'SMILES; round-tripping through SMILES will convert them to single bonds', LOST))

        for n, a in mol.atoms():
            if a.stereo is None:
                continue
            if n not in mol.stereogenic_tetrahedrons:
                continue  # allene or not a real stereocentre
            env = list(mol._bonds[n])
            s = mol._translate_tetrahedron_sign(n, env)
            pyramid = [mapping[x] for x in env]
            while len(pyramid) < 4:
                pyramid.append(-1)  # implicit hydrogen
            if s:
                pyramid[1], pyramid[2] = pyramid[2], pyramid[1]
            ig_mol.getAtom(mapping[n]).addStereocenter(
                _Indigo.ABS, pyramid[0], pyramid[1], pyramid[2], pyramid[3]
            )

    return ig_mol


def from_indigo(data, /, *, log=None):
    """
    Import an Indigo molecule as a chython 3 ``MoleculeContainer``.

    Aromaticity is taken verbatim from Indigo (order 4 goes in as order 4); chython does not re-run
    its own perception.  Tetrahedral stereo is translated from Indigo's pyramid.  Reportable losses:
    cis/trans (Indigo's bond-stereo values are derived from 2D coordinates and cannot be decoded back
    reliably), allene (no Indigo form), and OR/AND stereo groups, which arrive as plain parities with
    the group dropped.

    EVERY RECORD LANDS ON THE RETURNED MOLECULE'S `.log`, in stage `'interop'`, with nothing passed
    in; `log`, when given, receives a copy of the same records.

    :param data: an Indigo ``IndigoObject`` representing a molecule.
    :param log: optional list receiving a copy of the records put on the molecule.
    :raises UnconvertibleType: if `data` is not an Indigo IndigoObject.
    :raises ToolkitError: if Indigo itself raises an error during atom/bond iteration.
    """
    from indigo.indigo.indigo_object import IndigoObject as _IndigoObject
    from ..core import MoleculeContainer as _CoreMC
    from ..core._core import element_symbols as _element_symbols
    from ..exceptions import UnconvertibleType, ToolkitError

    if not isinstance(data, _IndigoObject):
        raise UnconvertibleType(
            f'from_indigo reads an Indigo IndigoObject, not {type(data).__qualname__!r}; '
            f'pass an Indigo molecule object or a chython container'
        )

    _syms = _element_symbols()
    _sym_to_elem = {s: i for i, s in enumerate(_syms) if s}

    mol = _CoreMC()
    idx_to_sid = {}
    records = []

    try:
        for a in data.iterateAtoms():
            idx = a.index()
            symbol = a.symbol()
            atomic = _sym_to_elem.get(symbol)
            if atomic is None:
                raise ToolkitError(
                    f'Indigo molecule contains unrecognised element symbol {symbol!r}'
                )
            charge = a.charge()
            isotope = a.isotope()
            radical = a.radical() != 0  # Indigo: 0=none, 101=singlet, 102=doublet, 103=triplet
            h = a.countImplicitHydrogens()
            sid = mol.add_atom(atomic, charge=charge, isotope=isotope, radical=radical, implicit_h=h)
            idx_to_sid[idx] = sid
    except (ToolkitError, UnconvertibleType):
        raise
    except Exception as exc:
        raise ToolkitError(f'Indigo raised during atom iteration: {exc}') from exc

    try:
        for b in data.iterateBonds():
            n = idx_to_sid[b.source().index()]
            m = idx_to_sid[b.destination().index()]
            order = b.bondOrder()
            mol.add_bond(n, m, order)
    except (ToolkitError,):
        raise
    except Exception as exc:
        raise ToolkitError(f'Indigo raised during bond iteration: {exc}') from exc

    from indigo import Indigo as _Indigo
    _EITHER = _Indigo.EITHER  # type 4: configuration not known

    for sc in data.iterateStereocenters():
        sc_type = sc.stereocenterType()
        if sc_type == _EITHER:
            continue  # unspecified configuration; leave parity=0

        pyramid = sc.stereocenterPyramid()
        anchor_idx = sc.index()
        anchor_sid = idx_to_sid[anchor_idx]

        # Indigo atom indices -> stable ids; -1 encodes implicit H (None in V3 refs)
        mapped = tuple(idx_to_sid[v] if v != -1 else None for v in pyramid)

        # Indigo reads its pyramid clockwise from the first entry, which is what the export above
        # writes for a stored parity of 1, so `mapped` read as a chython direction order *is* parity 1.
        # The probe asks the core which stored parity reads back that way instead of counting
        # inversions here, where a sign error would be a silently mirrored molecule.
        if not set_parity_by_probe(mol, anchor_sid, mapped, 1):
            # `stereo_units()` and not `stereogenic_units()`: an atom can be topologically equivalent
            # to a neighbour while every parity is still 0, and the post-automorphism subset would call
            # that atom "not a stereocentre".
            refs = next((u['refs'] for u in mol.stereo_units()
                         if u['anchor'] == anchor_sid and u['kind'] == SU_TETRA), None)
            if refs is None:
                records.append(LogRecord('indigo:tetrahedral-unit-not-found', (anchor_sid,),
                                         f'atom {anchor_sid}: Indigo reports a stereocenter but chython sees no '
                                         f'tetrahedral unit here; stereo not set', LOST))
            else:
                records.append(LogRecord('indigo:tetrahedral-pyramid-invalid', (anchor_sid,),
                                         f'atom {anchor_sid}: Indigo pyramid {pyramid!r} is not a permutation of '
                                         f'chython\'s directions {refs!r}; tetrahedral stereo not set', LOST))

        if sc_type in (2, 3):  # OR=2, AND=3 in Indigo
            group_name = 'OR' if sc_type == 2 else 'AND'
            records.append(LogRecord('indigo:stereo-group-converted-as-abs', (anchor_sid,),
                                     f'atom {anchor_sid}: Indigo stereo type {group_name} converted as ABS; '
                                     f'stereo group information is lost', LOST))

    has_ct_stereo = any(
        b.bondStereo() != 0
        for b in data.iterateBonds()
        if b.bondOrder() == 2
    )
    if has_ct_stereo:
        records.append(LogRecord('indigo:cis-trans-not-importable', (),
                                 'cis/trans stereo not imported: Indigo bond-stereo values (derived from 2D '
                                 'coordinates) cannot be reliably decoded to chython parity without coordinate '
                                 'geometry', LOST))

    records.append(LogRecord('indigo:coordinates-not-imported', (),
                             'coordinates are not imported; the source molecule\'s layout is dropped', LOST))
    deliver(mol, records, log)
    return mol
