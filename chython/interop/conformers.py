# -*- coding: utf-8 -*-
#
#  Copyright 2025, 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
3D conformer generation over RDKit's ETKDG and CDPKit's `ConformerGenerator`.

`generate_conformers` stores the generated geometry as models of the molecule and returns a count.  Each
engine adds its own hydrogens -- RDKit's `AddHs` appends them after the heavy atoms, CDPKit's
`prepareForConformerGeneration` completes them in place -- so the coordinates each engine hands back stay
keyed by the caller's atom numbers over the heavy atoms both engines keep first.
"""
from typing import Literal


def generate_conformers(mol, /, limit: int = 10, *, optimize: bool = False,
                        engine: Literal['rdkit', 'cdpkit'] | None = None,
                        **kwargs) -> int:
    """
    Generate 3D conformers for a molecule and store each as one of its models.

    Hydrogens are added by the engine so the geometry is generated for a complete structure, but only the
    heavy atoms of *mol* are placed.  The engines are ``'rdkit'`` (ETKDG) and ``'cdpkit'``
    (``ConformerGenerator``, https://pubs.acs.org/doi/10.1021/acs.jcim.3c00563).

    The generated set REPLACES any models the molecule carried, where XYZ frames append: a run is the
    whole answer, and two runs of ``limit=10`` would leave twenty models no field on a conformer could
    tell apart.  Every one of them states ``ext_index`` None -- nothing generated came from a file.
    ``SEG_XY`` is untouched: a generated geometry is not a layout.

    :param mol: a ``chython.core.MoleculeContainer``.
    :param limit: maximum number of conformers to generate.
    :param optimize: optimise with the MMFF94 force field (RDKit engine only).
    :param engine: override the engine chosen by ``chython.interop.config.conformer_engine``.
    :param kwargs: engine arguments.  CDPKit takes ``timeout`` (seconds, default 60), ``min_rmsd``
        (default .5) and ``energy_window`` (default 20); RDKit forwards everything to
        ``EmbedMultipleConfs``.
    :returns: how many conformers were stored, 0 when the engine produced nothing.
    :raises ValueError: if *engine* names an engine that is not implemented.
    """
    if engine is None:
        # read at call time and through `config`, so assigning the knob after import takes effect
        from . import config

        engine = config.conformer_engine

    if engine == 'rdkit':
        builder = _rdkit_conformers
    elif engine == 'cdpkit':
        builder = _cdpkit_conformers
    else:
        raise ValueError(f'no conformer generation engine named {engine!r}; use "rdkit" or "cdpkit"')

    # `atom_numbers` and not the engine's indices: both engines append their hydrogens after them
    coordinates = builder(mol, mol.atom_numbers, limit, optimize, kwargs)
    if not coordinates:
        # Nothing generated is nothing dropped: an engine that gave up is not a reason to lose the
        # models the caller already had.
        return 0

    existing = len(mol.conformers)
    if existing:
        # TWO SESSIONS, because one may not both drop and add: the drop shifts every index above it.
        with mol.edit():
            for i in range(existing):
                mol.drop_conformer(i)
    with mol.edit():
        for one in coordinates:
            model = mol.add_conformer()
            for n, (x, y, z) in one.items():
                mol.set_xyz(n, x, y, z, model=model)
    return len(coordinates)


def _rdkit_conformers(mol, heavy, limit, optimize, kwargs) -> list[dict]:
    """Embed with ETKDG.  RDKit keeps the exporter's atom order, so the heavy atoms lead."""
    from rdkit.Chem import AddHs
    from rdkit.Chem.AllChem import EmbedMultipleConfs, MMFFOptimizeMolecule

    from ._rdkit import to_rdkit

    # `AddHs` appends the hydrogens, so the heavy atoms keep their positions and the `zip` below can drop
    # the tail.  ETKDG without them embeds stereocentres with nothing to be chiral about.
    rmol = AddHs(to_rdkit(mol, keep_mapping=False, keep_hydrogens=False))
    ids = EmbedMultipleConfs(rmol, numConfs=limit, **kwargs)
    if optimize:
        for i in ids:
            MMFFOptimizeMolecule(rmol, confId=i)

    # `zip` stops at `heavy`, dropping the trailing hydrogen positions.
    return [{n: tuple(v) for n, v in zip(heavy, conf.GetPositions())}
            for conf in rmol.GetConformers() if conf.Is3D()]


def _cdpkit_conformers(mol, heavy, limit, optimize, kwargs) -> list[dict]:
    """
    Generate with CDPKit's ConformerGenerator.

    The molecule is built by `interop.cdpkit`, which sets the stereo descriptors directly; handing CDPKit
    an SDF proxy instead loses chirality, since without a 2D layout its wedges are ambiguous and the
    generator picks a handedness at random.  `optimize` is not offered: CDPKit's generator already
    returns force-field-minimised conformers.
    """
    from CDPL import Chem, ConfGen

    from ._cdpkit import to_cdpkit

    cmol = to_cdpkit(mol)
    # Index i is the i-th atom of `mol` (`interop.cdpkit`'s interface).  Taken before
    # `prepareForConformerGeneration`, which appends the hydrogens, so these indices stay valid.
    pos = {n: i for i, n in enumerate(mol.atom_numbers)}

    ConfGen.prepareForConformerGeneration(cmol)
    gen = ConfGen.ConformerGenerator()
    gen.settings.timeout = kwargs.get('timeout', 60) * 1000
    gen.settings.minRMSD = kwargs.get('min_rmsd', .5)
    gen.settings.energyWindow = kwargs.get('energy_window', 20.)
    gen.settings.maxNumOutputConformers = limit
    if gen.generate(cmol) != ConfGen.ReturnCode.SUCCESS:
        # a generator that gave up is an empty answer, not an exception
        return []

    gen.setConformers(cmol)
    atom_of = {n: cmol.getAtom(pos[n]) for n in heavy}
    return [{n: tuple(Chem.getConformer3DCoordinates(a, i)) for n, a in atom_of.items()}
            for i in range(gen.getNumConformers())]


__all__ = ['generate_conformers']
