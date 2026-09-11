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
CDPKit exporter (export only).  Reached through ``interop.cdpkit``.

CDPKit's InChI writer recomputes hydrogen counts from raw bond orders and treats order 4 as an
integer, so an aromatic molecule's formula comes back with hydrogens added (benzene -> C6H12); SMILES
is unaffected and is the oracle to use.  Atom order is part of the interface: 0-based ``cdp_mol.getAtom(i)`` is the
*i*-th atom of ``mol.atom_numbers``.
"""
from ..core import LogRecord, LOST, SU_CIS_TRANS, SU_TETRA
from ..exceptions import DirectionNotImplemented
from ._stereo import cis_trans_frame, kind_name


def _to_cdpkit_v3(mol, log, Chem):
    """Build a CDPKit BasicMolecule from a chython 3 MoleculeContainer."""
    from CDPL.Chem import AtomConfiguration, BondConfiguration, StereoDescriptor

    mol_cdp = Chem.BasicMolecule()
    idx_map = {}    # chython atom number → 0-based CDPKit index
    cdp_atoms = []  # CDPKit atom objects, indexed by position

    h_loss_logged = False

    for i, n in enumerate(mol.atom_numbers):
        a = mol_cdp.addAtom()
        anum = mol.element_of(n)
        Chem.setType(a, anum)
        Chem.setSymbol(a, Chem.getSymbolForType(a))

        charge = mol.charge_of(n)
        if charge:
            Chem.setFormalCharge(a, charge)

        isotope = mol.isotope_of(n)
        if isotope:
            Chem.setIsotope(a, isotope)

        if mol.radical_of(n):
            Chem.setRadicalType(a, Chem.RadicalType.DOUBLET)

        nH = mol.implicit_h_of(n)
        if nH is not None:
            Chem.setImplicitHydrogenCount(a, nH)
        elif not h_loss_logged:
            if log is not None:
                log.append(LogRecord('cdpkit:h-unknown-not-representable', (),
                                     'implicit hydrogen count is unknown for one or more atoms; '
                                     'CDPKit will infer hydrogen counts from valence', LOST))
            h_loss_logged = True

        idx_map[n] = i
        cdp_atoms.append(a)

    cdp_bonds = {}  # (n, m) both orderings → CDPKit bond object
    for n in mol.atom_numbers:
        for m in mol.neighbors_of(n):
            if m > n:
                order = mol.order_of(n, m)
                b = mol_cdp.addBond(idx_map[n], idx_map[m])
                Chem.setOrder(b, order)
                if order == 4:
                    Chem.setAromaticityFlag(b, True)
                    Chem.setAromaticityFlag(cdp_atoms[idx_map[n]], True)
                    Chem.setAromaticityFlag(cdp_atoms[idx_map[m]], True)
                cdp_bonds[n, m] = cdp_bonds[m, n] = b

    for unit in mol.stereo_units():
        kind = unit['kind']
        parity = unit['parity']
        if parity == 0:
            continue  # unset

        anchor = unit['anchor']
        refs = unit['refs']

        if kind == SU_TETRA:
            config = AtomConfiguration.R if parity == 1 else AtomConfiguration.S
            # a None ref is an implicit hydrogen, which CDPKit encodes as the central atom itself
            center_cdp = cdp_atoms[idx_map[anchor]]

            def _ref(r):
                return cdp_atoms[idx_map[r]] if r is not None else center_cdp

            sd = StereoDescriptor(config, _ref(refs[0]), _ref(refs[1]),
                                  _ref(refs[2]), _ref(refs[3]))
            Chem.setStereoDescriptor(center_cdp, sd)

        elif kind == SU_CIS_TRANS:
            # Framed through `_stereo`, and the parity is re-read through `translate_stereo` in the
            # frame actually written: refs positions may not be read directly (an unnamed slot may sit
            # inside a half) and a stored parity only means anything in the unit's own refs order.
            frame = cis_trans_frame(mol, anchor, refs)
            if frame is None:
                if log is not None:
                    log.append(LogRecord('cdpkit:cis-trans-frame-not-found', (anchor,),
                                         f'cis/trans stereo at atom {anchor}: cannot be framed on the graph; '
                                         f'skipped', LOST))
                continue
            bond = cdp_bonds.get((anchor, frame.partner))
            if bond is None:
                if log is not None:
                    log.append(LogRecord('cdpkit:cis-trans-bond-not-found', (anchor,),
                                         f'cis/trans stereo at atom {anchor}: the bond to its partner '
                                         f'{frame.partner} is not in the CDPKit molecule; skipped', LOST))
                continue
            try:
                framed = mol.translate_stereo(anchor, frame.order)
            except (KeyError, ValueError) as e:
                if log is not None:
                    log.append(LogRecord('cdpkit:cis-trans-config-not-readable', (anchor,),
                                         f'cis/trans stereo at atom {anchor}: chython will not read the '
                                         f'configuration in CDPKit\'s frame ({e}); skipped', LOST))
                continue
            config = BondConfiguration.TRANS if framed == 1 else BondConfiguration.CIS
            sd = StereoDescriptor(
                config,
                cdp_atoms[idx_map[frame.near]],
                cdp_atoms[idx_map[anchor]],
                cdp_atoms[idx_map[frame.partner]],
                cdp_atoms[idx_map[frame.far]],
            )
            Chem.setStereoDescriptor(bond, sd)

        elif log is not None:
            # every remaining kind, named from the core's own table so a new one cannot go unreported
            log.append(LogRecord('cdpkit:stereo-kind-not-representable', (anchor,),
                                 f'{kind_name(kind)} stereo at atom {anchor} cannot be represented in CDPKit; '
                                 f'skipped', LOST))

    Chem.calcBasicProperties(mol_cdp, False)
    return mol_cdp


def to_cdpkit(mol, /, *, log=None):
    """
    Export a chython container as a CDPKit ``BasicMolecule``.

    Atoms are emitted in the molecule's own iteration order.  Aromatic bonds (order 4) are stored with
    an aromatic flag; SMILES generation works, InChI generation does not — see the module note.
    Allene and atropisomer stereo cannot be represented and each skipped unit is logged.

    :param mol: a `chython.core.MoleculeContainer` to export.
    :param log: optional list; receives one string per reportable loss.
    :returns: a ``CDPL.Chem.BasicMolecule``.
    :raises: :class:`~chython.exceptions.UnconvertibleType` if *mol* is not a molecule container;
        :class:`~chython.exceptions.ToolkitError` if CDPKit raises unexpectedly.
    """
    from CDPL import Chem

    from ..core import MoleculeContainer as _V3Mol
    from ..exceptions import ToolkitError, UnconvertibleType

    if not isinstance(mol, _V3Mol):
        raise UnconvertibleType(f'{type(mol).__qualname__} has no CDPKit form: only a molecule '
                                f'container converts, and a query, a reaction and a CGR are not '
                                f'molecules')
    try:
        return _to_cdpkit_v3(mol, log, Chem)
    except Exception as exc:
        raise ToolkitError(str(exc)) from exc


def from_cdpkit(data, /, *, log=None):
    """
    Not built: chython reads no CDPKit molecule.  Raises ``DirectionNotImplemented``.
    """
    # When this direction is built: coordinates are not imported, so add a log line and extend
    # test_coordinate_honesty.py; and the records go to the returned molecule's `.log` through
    # `_records.deliver`, which test_log_delivery.py gets a row for.
    raise DirectionNotImplemented(
        'interop.cdpkit is export only; chython cannot read a CDPKit '
        'molecule. Write it to a file and read that, or use another toolkit'
    )
