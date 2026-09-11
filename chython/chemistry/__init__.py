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
"""Chemical knowledge as TSV in `tables/`, and the passes that apply it.

Imports `chython.core` and the standard library only, never the `chython` facade.  Importing it registers
`standardize()`, `canonicalize()` and friends on the core container by injection, `MoleculeContainer`
being a `cdef class`.  No pass runs on parse.

`perceive_bonds`, `saturate` and `expand_abbreviations` are the passes with no method.  The first two
are the two halves of building a molecule out of a coordinate file -- which pairs are bonded, then at
what order -- so their caller is whoever read that file, and neither runs on read.  The third reads what
a drawing wrote on an atom, which is a fact about a FILE and not about a structure, so it belongs beside
the reader that stored the alias rather than on every molecule.
"""
from ._abbreviations import expand_abbreviations
from ._canonicalize import canonicalize
from ._counts import (hydrogen_bond_acceptors_count, hydrogen_bond_donors_count,
                      rotatable_bonds_count)
from ._crippen import crippen_logp, crippen_mr
from ._hydrogens import explicify_hydrogens, implicify_hydrogens
from ._implicit import calc_implicit, check_valence
from ._isomers import standardize_isomers
from ._maccs import maccs_bit_set, maccs_keys
from ._perceive import perceive_bonds
from ._pharmacophore import pharmacophore_invariants
from ._protomers import neutralize
from ._qed import alert_count, qed, qed_properties
from ._residues import (RESIDUE_KINDS, ResidueTemplate, normalize_atom_name, residue_template,
                        residue_templates)
from ._resonance import fix_resonance
from ._salts import SaltComposition, decompose_salts, split_salts
from ._saturate import saturate
from ._smarts import SmartsSyntaxError, compile_smarts
from ._standardize import LogRecord, standardize
from ._tables import (ACID_ROLES, AbbreviationRow, AcidRow, Endpoint, RESONANCE_ROLES, Rule, SALT_ROLES,
                      SaltRow, abbreviation_row, abbreviations_rows, acids_rules, acids_rules_by_role,
                      acids_table_text, groups_rules,
                      metals_rules, read_table, resonance_rules, resonance_rules_by_role,
                      resonance_table_text, salts_rows, salts_rows_by_role, salts_species_keys,
                      salts_table_text, standardize_rules)
from ._tpsa import tpsa
from ..core._core import (_set_canonicalize_fn, _set_featurizer_fns, _set_hydrogens_fns,
                          _set_isomers_fn, _set_protomers_fn, _set_resonance_fn, _set_salts_fns,
                          _set_standardize_fn, _set_valence_fn)


__all__ = ['ACID_ROLES', 'AbbreviationRow', 'LogRecord', 'SALT_ROLES', 'SaltComposition',
           'abbreviation_row', 'abbreviations_rows', 'alert_count',
           'calc_implicit', 'canonicalize', 'check_valence', 'crippen_logp', 'crippen_mr',
           'decompose_salts', 'expand_abbreviations', 'explicify_hydrogens', 'fix_resonance',
           'hydrogen_bond_acceptors_count',
           'hydrogen_bond_donors_count', 'implicify_hydrogens', 'maccs_bit_set', 'maccs_keys',
           'neutralize', 'perceive_bonds', 'pharmacophore_invariants', 'qed', 'qed_properties',
           'rotatable_bonds_count',
           'saturate', 'split_salts', 'standardize', 'standardize_isomers', 'tpsa']

_set_standardize_fn(standardize)
_set_canonicalize_fn(canonicalize)
_set_hydrogens_fns(implicify_hydrogens, explicify_hydrogens)
_set_isomers_fn(standardize_isomers)
_set_valence_fn(check_valence)
_set_salts_fns(split_salts=split_salts, decompose_salts=decompose_salts)
_set_protomers_fn(neutralize)
_set_resonance_fn(fix_resonance)
_set_featurizer_fns(rotatable_bonds_count=rotatable_bonds_count,
                    hydrogen_bond_donors_count=hydrogen_bond_donors_count,
                    hydrogen_bond_acceptors_count=hydrogen_bond_acceptors_count,
                    tpsa=tpsa, crippen_logp=crippen_logp, crippen_mr=crippen_mr, qed=qed,
                    maccs_keys=maccs_keys, maccs_bit_set=maccs_bit_set,
                    pharmacophore_invariants=pharmacophore_invariants)
