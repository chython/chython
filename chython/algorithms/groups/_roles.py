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
from lazy_object_proxy import Proxy


# Role glossary: role name -> [(functional-group pattern name, cap template), ...].
# The cap template is a reactor product SMARTS. `[A:N]` preserves the source
# atom; the pattern's `:100` (and orphaned atoms) are deleted; the attachment
# atom is capped with a new `[At:20]`.
#
# A role groups functional groups that expose the same coupling handle, split by
# the hybridization of the attachment carbon (aryl / alkyl / alkenyl-sp2 /
# alkynyl-sp) so complementary roles pair cleanly into coupling enumerators.
# Some functional groups appear under several roles on purpose (e.g. a phenol is
# both an O-nucleophile keeping its oxygen for Mitsunobu and a deoxygenative
# C-coupling handle); downstream dedup keeps the roles distinct.
#
# Cap-template notes tied to the shape of the reused `_functional.py` patterns:
#   * alkenyl patterns map the second sp2 carbon as `:2` -> keep it: `[A:1](=[A:2])...`
#   * alkynyl halide/boron patterns do NOT map the far sp carbon -> `[A:1]-[At]`
#     (the triple bond is preserved from the source molecule); alkynyl_silane
#     DOES map `:2`, so it keeps `[A:1](#[A:2])...`.
#   * geometry carboxylic-acid patterns map the R-carbon as `:3` -> acyl caps
#     keep it (`[A:3]-[A:1](=[A:2])...`) or the R side is orphan-deleted.
#
# This glossary is maintainer-curated chemistry: extend/edit `roles` here without
# touching the consuming code.
roles = {
    # --- electrophiles: halide leaving groups (Cl/Br/I), fluoride separate ---
    'aryl_halide': [
        ('aryl_chloride', '[A:1]-[At:20]'),
        ('aryl_bromide',  '[A:1]-[At:20]'),
        ('aryl_iodide',   '[A:1]-[At:20]'),
    ],
    'alkyl_halide': [
        ('alkyl_chloride', '[A:1]-[At:20]'),
        ('alkyl_bromide',  '[A:1]-[At:20]'),
        ('alkyl_iodide',   '[A:1]-[At:20]'),
    ],
    'alkenyl_halide': [
        ('alkenyl_chloride', '[A:1](=[A:2])-[At:20]'),
        ('alkenyl_bromide',  '[A:1](=[A:2])-[At:20]'),
        ('alkenyl_iodide',   '[A:1](=[A:2])-[At:20]'),
    ],
    'alkynyl_halide': [
        ('alkynyl_chloride', '[A:1]-[At:20]'),
        ('alkynyl_bromide',  '[A:1]-[At:20]'),
        ('alkynyl_iodide',   '[A:1]-[At:20]'),
    ],
    # fluoride kept separate from the other halides
    'aryl_fluoride':    [('aryl_fluoride',    '[A:1]-[At:20]')],
    'alkyl_fluoride':   [('alkyl_fluoride',   '[A:1]-[At:20]')],
    'alkenyl_fluoride': [('alkenyl_fluoride', '[A:1](=[A:2])-[At:20]')],
    'alkynyl_fluoride': [('alkynyl_fluoride', '[A:1]-[At:20]')],
    # sulfonate pseudohalide leaving groups
    'aryl_sulfonate': [
        ('aryl_triflate',  '[A:1]-[At:20]'),
        ('aryl_mesylate',  '[A:1]-[At:20]'),
        ('aryl_tosylate',  '[A:1]-[At:20]'),
    ],
    'alkyl_sulfonate': [
        ('alkyl_triflate', '[A:1]-[At:20]'),
        ('alkyl_mesylate', '[A:1]-[At:20]'),
        ('alkyl_tosylate', '[A:1]-[At:20]'),
    ],

    # --- nucleophilic carbon: boron + organometal transmetalation donors ---
    'aryl_boron': [
        ('aryl_boronic_acid',   '[A:1]-[At:20]'),
        ('aryl_boronic_ester',  '[A:1]-[At:20]'),
        ('aryl_molander_salt',  '[A:1]-[At:20]'),
    ],
    'alkyl_boron': [
        ('alkyl_boronic_acid',   '[A:1]-[At:20]'),
        ('alkyl_boronic_ester',  '[A:1]-[At:20]'),
        ('alkyl_molander_salt',  '[A:1]-[At:20]'),
    ],
    'alkenyl_boron': [
        ('alkenyl_boronic_acid',  '[A:1](=[A:2])-[At:20]'),
        ('alkenyl_boronic_ester', '[A:1](=[A:2])-[At:20]'),
        ('alkenyl_molander_salt', '[A:1](=[A:2])-[At:20]'),
    ],
    'alkynyl_boron': [
        ('alkynyl_boronic_acid',  '[A:1]-[At:20]'),
        ('alkynyl_boronic_ester', '[A:1]-[At:20]'),
        ('alkynyl_molander_salt', '[A:1]-[At:20]'),
    ],
    'aryl_magnesium':    [('aryl_grignard',    '[A:1]-[At:20]')],
    'alkyl_magnesium':   [('alkyl_grignard',   '[A:1]-[At:20]')],
    'alkenyl_magnesium': [('alkenyl_grignard', '[A:1](=[A:2])-[At:20]')],
    'aryl_zinc':         [('aryl_zinc',        '[A:1]-[At:20]')],
    'alkyl_zinc':        [('alkyl_zinc',       '[A:1]-[At:20]')],
    'alkenyl_zinc':      [('alkenyl_zinc',     '[A:1](=[A:2])-[At:20]')],
    'aryl_stannane':     [('aryl_stannane',    '[A:1]-[At:20]')],
    'alkyl_stannane':    [('alkyl_stannane',   '[A:1]-[At:20]')],
    'alkenyl_stannane':  [('alkenyl_stannane', '[A:1](=[A:2])-[At:20]')],
    'aryl_silane':       [('aryl_silane',      '[A:1]-[At:20]')],
    'alkenyl_silane':    [('alkenyl_silane',   '[A:1](=[A:2])-[At:20]')],
    'alkynyl_silane':    [('alkynyl_silane',   '[A:1](#[A:2])-[At:20]')],

    # --- acyl handle: keep the carbonyl O and the R-carbon (:3) ---
    'alkyl_acyl':    [('alkyl_carboxylic_acid',    '[A:3]-[A:1](=[A:2])-[At:20]')],
    'aryl_acyl':     [('aryl_carboxylic_acid',     '[A:3]-[A:1](=[A:2])-[At:20]')],
    # alkenyl acid pattern maps the beta-vinyl carbon as :4 -> keep it too
    'alkenyl_acyl':  [('alkenyl_carboxylic_acid',  '[A:3](=[A:4])-[A:1](=[A:2])-[At:20]')],
    'alkynyl_acyl':  [('alkynyl_carboxylic_acid',  '[A:3]-[A:1](=[A:2])-[At:20]')],
    'acyl_halide': [
        ('acyl_chloride', '[A:1](=[A:2])-[At:20]'),
        ('acyl_bromide',  '[A:1](=[A:2])-[At:20]'),
        ('acyl_fluoride', '[A:1](=[A:2])-[At:20]'),
    ],
    # carbamoyl electrophile (R2N-C(=O)-X -> urea/carbamate): cap the carbonyl carbon,
    # keep the carbonyl O (:2) and the amine N (:3); the halide leaving group is deleted.
    'carbamoyl_halide': [
        ('carbamoyl_chloride', '[A:3]-[A:1](=[A:2])-[At:20]'),
        ('carbamoyl_fluoride', '[A:3]-[A:1](=[A:2])-[At:20]'),
    ],

    # --- heteroatom nucleophiles: cap the heteroatom, keep N/O/S ---
    'alkyl_amine': [
        ('primary_amine',   '[A:1](-[A:2])-[At:20]'),
        ('secondary_amine', '[A:1](-[A:2])(-[A:3])-[At:20]'),
    ],
    'aryl_amine': [
        ('primary_aniline',   '[A:1](-[A:2])-[At:20]'),
        ('secondary_aniline', '[A:1](-[A:2])(-[A:3])-[At:20]'),
    ],
    # amide N-H as Goldberg/Buchwald N-nucleophile: keep N(:1) + carbonyl C(:2)=O(:3)
    'amide_nitrogen': [
        ('primary_amide',   '[A:1](-[A:2]=[A:3])-[At:20]'),
        ('secondary_amide', '[A:1](-[A:2]=[A:3])-[At:20]'),
    ],
    # amidine NH2 as Buchwald-Hartwig N-nucleophile: cap the primary N(:1), keep the
    # amidine carbon (:2); the far imine N survives as an orphan hanging off :2.
    'amidine_nitrogen': [('primary_amidine_amine', '[A:1](-[A:2])-[At:20]')],
    # azole N-H as N-arylation nucleophile: mirror _reactions.py ullmann_pyrrole --
    # the handle attaches to the h0 (pyridine-type) nitrogen (pyrazole :2, imidazole
    # :3), NOT the matched h1 NH (:1). Redraw ring bonds aromatic (:) so the ring
    # survives; the stale h1 on the far nitrogen heals via ignore_pyrrole_hydrogen.
    'azole_nitrogen': [
        ('pyrrole',   '[A:1]-[At:20]'),
        ('pyrazole',  '[A:1]:[A:2]-[At:20]'),
        ('imidazole', '[A:1]:[A:2]:[A:3]-[At:20]'),
    ],
    'alkyl_thiol': [('thiol', '[A:1](-[A:2])-[At:20]')],
    'aryl_thiol':  [('aryl_thiol', '[A:1](-[A:2])-[At:20]')],
    # O-nucleophiles: keep the oxygen (-> R-O[At] / Ar-O[At] / R-C(=O)-O[At])
    'alkyl_hydroxyl': [
        ('primary_alcohol',   '[A:1](-[A:2])-[At:20]'),
        ('secondary_alcohol', '[A:1](-[A:2])-[At:20]'),
        ('tertiary_alcohol',  '[A:1](-[A:2])-[At:20]'),
    ],
    # phenol keeps its oxygen directly; pyridone (azinone C:1=O:2, ring N:3) exposes
    # the same aryl-O handle via its hydroxyazine tautomer -- the ring aromatizes
    # (C=O -> C-O[At], C-N -> C=N), mirroring the mitsunobu/ullmann_phenol O-couplings.
    'aryl_hydroxyl': [
        ('phenol',  '[A:1](-[A:2])-[At:20]'),
        ('azinone', '[A:1](=[A:3])-[A:2]-[At:20]'),
    ],
    'acid_hydroxyl': [('carboxylic_acid', '[A:1](=[A:2])-[A:100]-[At:20]')],

    # terminal alkyne C-H as Sonogashira nucleophile: cap the terminal sp carbon,
    # keep the whole triple bond (:1#:2).
    'alkynyl_terminal': [('terminal_alkyne', '[At:20]-[A:1]#[A:2]')],

    # sulfonyl electrophile (RSO2X -> sulfonamide/sulfonate ester): cap the sulfur,
    # keep both S=O (:2, :3); the halide leaving group (:100) is orphan-deleted.
    'sulfonyl': [
        ('sulfonyl_chloride', '[A:1](=[A:2])(=[A:3])-[At:20]'),
        ('sulfonyl_fluoride', '[A:1](=[A:2])(=[A:3])-[At:20]'),
    ],

    # --- deoxygenative coupling: drop the O, cap the carbon ---
    'alkyl_deoxy': [
        ('primary_alcohol',   '[A:2]-[At:20]'),
        ('secondary_alcohol', '[A:2]-[At:20]'),
        ('tertiary_alcohol',  '[A:2]-[At:20]'),
    ],
    'aryl_deoxy': [('phenol', '[A:2]-[At:20]')],

    # --- reductive amination C-side: drop the carbonyl O, cap the carbon (:1) ---
    'carbonyl_electrophile': [
        ('aldehyde', '[A:1]-[At:20]'),
        ('ketone',   '[A:1]-[At:20]'),
    ],

    # --- decarboxylative coupling: drop the whole COOH, cap the R-carbon (:3) ---
    'alkyl_decarboxy':   [('alkyl_carboxylic_acid',    '[A:3]-[At:20]')],
    'aryl_decarboxy':    [('aryl_carboxylic_acid',     '[A:3]-[At:20]')],
    'alkenyl_decarboxy': [('alkenyl_carboxylic_acid',  '[A:3](=[A:4])-[At:20]')],
    'alkynyl_decarboxy': [('alkynyl_carboxylic_acid',  '[A:3]-[At:20]')],

    # --- deaminative coupling: drop the N, cap the alpha carbon (:2) ---
    # sp3 amines couple via Katritzky pyridinium salts; anilines via diazonium.
    'alkyl_deamino': [('primary_amine',   '[A:2]-[At:20]')],
    'aryl_deamino':  [('primary_aniline', '[A:2]-[At:20]')],
}


def _transformers():
    from ._functional import rules as functional_rules
    from ...reactor import Transformer
    from ... import smarts

    compiled = {}
    for role_name, entries in roles.items():
        built = []
        for fg_name, cap_template in entries:
            t = Transformer(functional_rules[fg_name], smarts(cap_template),
                            delete_atoms=True, automorphism_filter=True, canonicalize=True,
                            ignore_pyrrole_hydrogen=True)
            built.append((fg_name, t))
        compiled[role_name] = built
    return compiled


transformers = Proxy(_transformers)


__all__ = ['roles', 'transformers']
