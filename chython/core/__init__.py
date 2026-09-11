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
from pathlib import Path
from ._core import (Atom, AutomorphismBudgetExceeded, Bond, MoleculeContainer, QueryContainer,
                    WEDGE_DOWN, WEDGE_EITHER, WEDGE_NONE, WEDGE_UP,
                    STEREO_ABS, STEREO_AND, STEREO_OR, STEREO_UNSPECIFIED,
                    # the stereo-unit kinds, because `stereo_units()` reports `kind` as an
                    # integer and a consumer switching on it should not have to spell the
                    # integer.  SU_HELICAL is reserved and never produced; it is here so that
                    # "not this one" can be written as a name too.
                    SU_TETRA, SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER, SU_HELICAL,
                    H_IMPLICIT_MAX, H_UNKNOWN, R_INDEX_MAX,
                    ich_load_library, inchi_library_loaded, inchi, inchikey,
                    molecule_to_inchi, molecule_to_inchikey, inchi_to_molecule,
                    _ich_set_kekule_fn,
                    isotope_data, isotope_offsets_table, isotope_counts_table,
                    # the SMILES writer.  `write_smiles` is what `str(mol)`, `format(mol, spec)`
                    # and `mol.smiles` all call, exported because a caller writing a million
                    # records wants the function without an attribute lookup per record, and
                    # `normalize_smiles_spec` because a caller CACHING those strings needs the
                    # spec resolved to the behaviour it selects rather than to its spelling.
                    write_smiles, normalize_smiles_spec, detached_smiles, DetachedSmiles,
                    # and the reaction writer, which is what `str(rxn)`, `format(rxn, spec)` and
                    # `rxn.smiles` all call.  It is a function and not a method for the reason
                    # `reaction.py` gives: aggregating one CXSMILES tail across three sides is the
                    # writer's machinery, and the container holds no formatting of its own.
                    write_reaction_smiles,
                    # `sticky_smiles` is the glue-two-strings-together fragment, exported
                    # because consumers outside this repository call the method it backs.
                    sticky_smiles,
                    # and the reader, which has the same reason plus one of its own: it is the
                    # only entry point that takes a `log`, so a pipeline that needs to see what
                    # was repaired cannot go through `smiles()` or any method.
                    read_smiles,
                    # and the reaction reader over it.  `read_smiles` returns a reaction when the
                    # string has an arrow, so this is not the only door -- it is the STRICT one,
                    # for the caller whose corpus is reactions and who wants the wrong shape to
                    # fail rather than to come back as a molecule.
                    read_reaction_smiles,
                    # and the SMARTS reader beside it, for the same reason and one more: it is
                    # the whole of the query front end, so a caller building patterns has no
                    # other entry point to reach for
                    read_smarts, IncorrectSmarts,
                    # and the SMIRKS reader over it.  `read_smirks` is the ONLY way to build a
                    # `ReactionTemplate` -- the class is exported for `isinstance` and for the
                    # type name in a traceback, and its `__init__` refuses -- because a template
                    # is a notation, and a second construction path is a second notation.
                    read_smirks, ReactionTemplate, IncorrectSmirks,
                    # the ML layout object: built once per dataset, reused across every view call
                    TensorEncoding,
                    # the per-atom view: element, hydrogens, heavy degree, distances as int32 arrays
                    StateView, mol_state_view,
                    # the per-atom before/after view over the union graph; a molecule is before==after
                    TransitionView, mol_transition_view, reaction_transition_view,
                    # the two representation changes, as functions beside the methods
                    kekule, KekuleResult, thiele, ThieleResult,
                    # the legacy pach codec.  `pach_load` is exported and `MoleculeContainer.unpack`
                    # is not enough on its own, because the method has to return a molecule or
                    # raise, and a caller walking a store of forty thousand records needs the door
                    # that reports a damaged record instead of ending the loop.
                    pach_load, pach_dump, pach_record_length)
# The reaction container, which is Python and not part of the extension: it holds three tuples and a
# title, so there is no loop in it for C to make faster.  See its module docstring, and note that the
# import comes AFTER `._core` because it imports from there.
from .reaction import (MappingResult, ReactionContainer, ReactionModelingView,
                       # the reaction-level pach codec, beside the molecule-level one above and for
                       # the same reason: `reaction_pach_load` reports a damaged record where
                       # `ReactionContainer.unpack` has to raise, and a loop over a store of packed
                       # reactions needs the door that does not end the loop.
                       reaction_pach_dump, reaction_pach_load)
# The bidirectional short doors, LAST because they import from both of the above: `smiles` is the
# spelling the whole tree uses for the SMILES reader and writer alike, and `pach`/`unpach`/`unpack`
# the same for the wire format.  Beside them and not instead of them, because a direction-stating
# function is what a loop calls: `read_smiles` takes the `log`, `pach_load` reports instead of
# raising, and neither is reachable through a door that has to decide which one was meant.
from ._facade import pach, smiles, unpach, unpack
# `LogRecord` is here and not in `chemistry` because it is not only `chemistry`'s: the SMIRKS patcher
# reports the parities it dropped in the same shape and lives in the extension, which cannot import
# upwards.  `chython.chemistry` re-exports the name, so both spellings are one class.
from ._log import INFO, LOST, REFUSED, REPAIRED, Log, LogRecord, recording


__all__ = ['Atom', 'AutomorphismBudgetExceeded', 'Bond', 'Log', 'LogRecord', 'recording',
           'INFO', 'REPAIRED', 'LOST', 'REFUSED',
           'MoleculeContainer', 'QueryContainer',
           'MappingResult', 'ReactionContainer', 'ReactionModelingView',
           'WEDGE_DOWN', 'WEDGE_EITHER', 'WEDGE_NONE', 'WEDGE_UP',
           'STEREO_ABS', 'STEREO_AND', 'STEREO_OR', 'STEREO_UNSPECIFIED',
           'SU_TETRA', 'SU_CIS_TRANS', 'SU_ALLENE', 'SU_ATROPISOMER', 'SU_HELICAL',
           'H_IMPLICIT_MAX', 'H_UNKNOWN', 'R_INDEX_MAX',
           'ich_load_library', 'inchi_library_loaded', 'inchi', 'inchikey',
           'molecule_to_inchi', 'molecule_to_inchikey', 'inchi_to_molecule',
           '_ich_set_kekule_fn',
           'isotope_data', 'isotope_offsets_table', 'isotope_counts_table',
           'write_smiles', 'write_reaction_smiles', 'normalize_smiles_spec', 'detached_smiles',
           'DetachedSmiles',
           'sticky_smiles', 'read_smiles', 'read_reaction_smiles', 'read_smarts', 'IncorrectSmarts',
           'read_smirks', 'ReactionTemplate', 'IncorrectSmirks',
           'TensorEncoding', 'StateView', 'mol_state_view',
           'TransitionView', 'mol_transition_view', 'reaction_transition_view',
           'kekule', 'KekuleResult', 'thiele', 'ThieleResult',
           'pach_load', 'pach_dump', 'pach_record_length',
           'reaction_pach_dump', 'reaction_pach_load',
           'smiles', 'pach', 'unpach', 'unpack']


def _auto_load_libinchi():
    """Try to load the bundled libinchi at import time.

    Looks for the library inside this package: chython/core/libinchi.so (Linux), libinchi.dylib
    (macOS), libinchi.dll (Windows).  Falls back silently if not found so that the module still
    imports; calling molecule_to_inchi() or inchi_to_molecule() will then raise ImportError with a
    clear message.

    The binary lives in `core/` and not beside a Python wrapper, because `core` is its only
    consumer: the bridge is `_inchi.pxi`, compiled into this package's one extension.
    """
    import sys
    _base = Path(__file__).parent
    _candidates = {
        'darwin': [_base / 'libinchi.dylib'],
        'win32':  [_base / 'libinchi.dll',  _base / 'libinchi.so'],
        'linux':  [_base / 'libinchi.so',   _base / 'libinchi.so.1'],
    }
    _platform = sys.platform
    for _prefix, _names in _candidates.items():
        if _platform.startswith(_prefix):
            for _p in _names:
                if _p.exists():
                    ich_load_library(str(_p))
                    return
            break


_auto_load_libinchi()
