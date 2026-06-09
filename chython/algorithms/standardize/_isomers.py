# -*- coding: utf-8 -*-
#
#  Copyright 2021-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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


def _rules():
    """
    Rules for canonical positioning of charges and hydrogens.

    Rule format: (query, rule_type)

    rule_type:
      'charge_fixed'         - deterministic charge move from atom_1 to atom_2
      'charge_fixed_bridge'  - deterministic charge move from atom_3 (bridge) to atom_2
      'charge_morgan'        - Morgan-dependent charge (atom_1 currently charged)
      'charge_morgan_bridge' - Morgan-dependent charge (atom_3 bridge currently charged)
      'ferrocene'            - ferrocene Cp ring charge
      'tautomer_fixed'       - deterministic H move
      'tautomer_morgan'      - Morgan-dependent H (atom_1 is donor)
      'tautomer_morgan_donor' - Morgan-dependent H (atom_3 is donor)
      'amidine'              - global H+bond resolution

    Charge rules encode atom validity constraints in SMARTS:
      - D2;h1: pyrrole-like NH+ (2 aromatic bonds + 1 implicit H)
      - D3([A]): substituted N+ candidate with any single-bonded atom

    Each charge rule has 4 variants for atom_1/atom_2 combinations:
      HH: both D2;h1
      HA: atom_1 D2;h1, atom_2 D3
      AH: atom_1 D3, atom_2 D2;h1
      AA: both D3
    """
    from ... import smarts

    rules = []

    # === charge_fixed: deterministic charge positioning ===
    # These rules match fused 5+5 and 5+6 ring systems where the canonical charge position is known.
    # charge_fixed: atom_1 is the current charge holder; charge moves to atom_2.
    # charge_fixed_bridge: atom_3 (bridge N) is the current charge holder; charge moves to atom_2.

    # --- Rule 1: N1C=CN2[NH+]=CC=C12 -> N1C=CC2=[NH+]C=CN12 ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[N;r5:3]:2:[C;r5](:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:2):[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[N;r5:3]:2:[C;r5](:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:2):[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[N;r5:3]:2:[C;r5](:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:2):[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[N;r5:3]:2:[C;r5](:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:2):[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))

    # --- Rule 2: N1C=C[N+]2=C1C=CN2 -> N1C=CC2=[NH+]C=CN12 ---
    # atom_3 is the charged bridge N (unconstrained D - can be D3 all-aromatic)
    # HH variant
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:[C;r5]:1:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # HA variant
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:[C;r5]:1:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # AH variant
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:[C;r5]:1:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # AA variant
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:[C;r5]:1:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))

    # --- Rule 3: N1C=C2C=CN[N+]2=C1 -> N1C=CC2=C[NH+]=CN12 ---
    # atom_3 is the charged bridge N (unconstrained D)
    # HH variant
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:[C;r5]:1:[C;r5]:[N;r5;D2;h1:2]:[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # HA variant
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:[C;r5]:1:[C;r5]:[N;r5;D3:2]([A]):[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # AH variant
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:[C;r5]:1:[C;r5]:[N;r5;D2;h1:2]:[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # AA variant
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:[C;r5]:1:[C;r5]:[N;r5;D3:2]([A]):[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))

    # --- Rule 4: N1C=C2NC=C[N+]2=C1 -> N1C=CN2C=[NH+]C=C12 ---
    # atom_3 is the charged bridge N (unconstrained D)
    # HH variant
    q = smarts('[N;a;r5;+:3]:1:2:[C;r5](:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:1):[C;r5]:[N;r5;D2;h1:2]:[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # HA variant
    q = smarts('[N;a;r5;+:3]:1:2:[C;r5](:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:1):[C;r5]:[N;r5;D3:2]([A]):[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # AH variant
    q = smarts('[N;a;r5;+:3]:1:2:[C;r5](:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:1):[C;r5]:[N;r5;D2;h1:2]:[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))
    # AA variant
    q = smarts('[N;a;r5;+:3]:1:2:[C;r5](:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:1):[C;r5]:[N;r5;D3:2]([A]):[C;r5]:2')
    rules.append((q, 'charge_fixed_bridge'))

    # --- Rule 5: N1C2=[NH+]C=CC2=CC=C1 -> N1C=CC2=CC=C[NH+]=C12 ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r6]:2:[N;r6;D2;h1:2]:[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r6]:2:[N;r6;D3:2]([A]):[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r6]:2:[N;r6;D2;h1:2]:[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r6]:2:[N;r6;D3:2]([A]):[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))

    # --- Rule 6: N1C=C2C(=CC=[NH+]2)C=C1 -> N1C=CC2=CC=[NH+]C=C12 ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r6]:2:[C;r6]:[N;r6;D2;h1:2]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r6]:2:[C;r6]:[N;r6;D3:2]([A]):[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r6]:2:[C;r6]:[N;r6;D2;h1:2]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r6]:2:[C;r6]:[N;r6;D3:2]([A]):[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))

    # --- Rule 7: C1=CC=2C(=[NH+]1)C=CNC=2 -> N1C=CC2=C[NH+]=CC=C12 ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r6]:2:[C;r6]:[C;r6]:[N;r6;D2;h1:2]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r6]:2:[C;r6]:[C;r6]:[N;r6;D3:2]([A]):[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r6]:2:[C;r6]:[C;r6]:[N;r6;D2;h1:2]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r6]:2:[C;r6]:[C;r6]:[N;r6;D3:2]([A]):[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))

    # --- Rule 8: C=1C=[NH+]C=2C=1NC=CC=2 -> N1C=CC2=[NH+]C=CC=C12 ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r6]:2:[C;r6]:[C;r6]:[C;r6]:[N;r6;D2;h1:2]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r6]:2:[C;r6]:[C;r6]:[C;r6]:[N;r6;D3:2]([A]):[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r6]:2:[C;r6]:[C;r6]:[C;r6]:[N;r6;D2;h1:2]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r6]:2:[C;r6]:[C;r6]:[C;r6]:[N;r6;D3:2]([A]):[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_fixed'))

    # --- Rule 9: C1=2C=[NH+]C=C1NC=CC=2 -> N1C=C2C=CC=[NH+]C2=C1 ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r5]:[C;r6]:2:[N;r6;D2;h1:2]:[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r5]:[C;r6]:2:[N;r6;D3:2]([A]):[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r5]:[C;r6]:2:[N;r6;D2;h1:2]:[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r5]:[C;r6]:2:[N;r6;D3:2]([A]):[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, 'charge_fixed'))

    # --- Rule 10: C1=2C=CNC=C1C=[NH+]C=2 -> N1C=C2C=C[NH+]=CC2=C1 ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r5]:[C;r6]:2:[C;r6]:[N;r6;D2;h1:2]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r5]:[C;r6]:2:[C;r6]:[N;r6;D3:2]([A]):[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r5]:[C;r6]:2:[C;r6]:[N;r6;D2;h1:2]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, 'charge_fixed'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r5]:[C;r6]:2:[C;r6]:[N;r6;D3:2]([A]):[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, 'charge_fixed'))

    # === tautomer_fixed: deterministic H positioning ===

    # 1,2,4-triazole: N1C=NC=N1 -> N1C=NN=C1
    q = smarts('[N;a;r5;D2;h1;x1:2]:1:[N;r5;D2;h0;x1]:[C;r5]:[N;r5;D2;h0;x0:1]:[C;r5]:1')
    rules.append((q, 'tautomer_fixed'))

    # tetrazole: N1N=CN=N1 -> N1C=NN=N1
    q = smarts('[N;a;r5;D2;h0;x1:1]:1:[N;r5;D2;h1;x2:2]:[N;r5;D2;h0;x2]:[N;r5;D2;h0;x1]:[C;r5]:1')
    rules.append((q, 'tautomer_fixed'))

    # === charge_morgan: Morgan-dependent charge positioning ===
    # Morgan ordering decides which atom gets the charge.

    # --- Rule M1: indazole+ fused (5+5 ring, N-N bridge) ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[N;D3]:2:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:[C;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[N;D3]:2:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:[C;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[N;D3]:2:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:[C;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[N;D3]:2:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:[C;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))

    # --- Rule M2: indazole+ fused (bridge variant, charge on bridge N) ---
    # HH variant for atom_1/atom_2
    q = smarts('[N;a;D3;r5;+:3]:1:2:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:[C;D3]:1:[C;r5]:[C;r5]:[N;r5;D2;h1:2]:2')
    rules.append((q, 'charge_morgan_bridge'))
    # HA variant
    q = smarts('[N;a;D3;r5;+:3]:1:2:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:[C;D3]:1:[C;r5]:[C;r5]:[N;r5;D3:2]([A]):2')
    rules.append((q, 'charge_morgan_bridge'))
    # AH variant
    q = smarts('[N;a;D3;r5;+:3]:1:2:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:[C;D3]:1:[C;r5]:[C;r5]:[N;r5;D2;h1:2]:2')
    rules.append((q, 'charge_morgan_bridge'))
    # AA variant
    q = smarts('[N;a;D3;r5;+:3]:1:2:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:[C;D3]:1:[C;r5]:[C;r5]:[N;r5;D3:2]([A]):2')
    rules.append((q, 'charge_morgan_bridge'))

    # --- Rule M3: imidazole fused (5+5, C-N-C bridge) ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;D3]:2:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:[N;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;D3]:2:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:[N;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;D3]:2:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:[N;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;D3]:2:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:[N;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))

    # --- Rule M4: imidazole fused (bridge variant, charge on bridge) ---
    # HH variant
    q = smarts('[N;a;D3;r5;+:3]:1:2:[C;D3](:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:1):[N;r5;D2;h1:2]:[C;r5]:[C;r5]:2')
    rules.append((q, 'charge_morgan_bridge'))
    # HA variant
    q = smarts('[N;a;D3;r5;+:3]:1:2:[C;D3](:[N;r5;D2;h1:1]:[C;r5]:[C;r5]:1):[N;r5;D3:2]([A]):[C;r5]:[C;r5]:2')
    rules.append((q, 'charge_morgan_bridge'))
    # AH variant
    q = smarts('[N;a;D3;r5;+:3]:1:2:[C;D3](:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:1):[N;r5;D2;h1:2]:[C;r5]:[C;r5]:2')
    rules.append((q, 'charge_morgan_bridge'))
    # AA variant
    q = smarts('[N;a;D3;r5;+:3]:1:2:[C;D3](:[N;r5;D3:1]([A]):[C;r5]:[C;r5]:1):[N;r5;D3:2]([A]):[C;r5]:[C;r5]:2')
    rules.append((q, 'charge_morgan_bridge'))

    # --- Rule M5: imidazolium (simple 5-ring) ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r5]:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[C;r5]:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r5]:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[C;r5]:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))

    # --- Rule M6: pyrazolium (simple 5-ring) ---
    # HH variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # HA variant
    q = smarts('[N;a;r5;D2;h1;+:1]:1:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # AH variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[N;r5;D2;h1:2]:[C;r5]:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))
    # AA variant
    q = smarts('[N;a;r5;D3;+:1]([A]):1:[N;r5;D3:2]([A]):[C;r5]:[C;r5]:[C;r5]:1')
    rules.append((q, 'charge_morgan'))

    # === ferrocene ===
    q = smarts('[C;a;r5;-:1]:1:[C;a;r5:2]:[C;a;r5:3]:[C;a;r5:4]:[C;a;r5:5]:1')
    rules.append((q, 'ferrocene'))

    # === tautomer_morgan: Morgan-dependent H positioning ===

    # pyrazole
    q = smarts('[N;a;r5;D2;h1;x1:1]:1:[N;r5;D2;h0;x1:2]:[C;r5]:[C;r5]:[C;r5]:1')
    rules.append((q, 'tautomer_morgan'))

    # imidazole
    q = smarts('[N;a;r5;D2;h1;x0:1]:1:[C;r5]:[N;r5;D2;h0;x0:2]:[C;r5]:[C;r5]:1')
    rules.append((q, 'tautomer_morgan'))

    # 1,2,3-triazole
    q = smarts('[N;a;r5;D2;h1;x1:1]:1:[N;r5;D2;h0;x2]:[N;r5;D2;h0;x1:2]:[C;r5]:[C;r5]:1')
    rules.append((q, 'tautomer_morgan'))

    # 1,2,3-triazole (donor variant)
    q = smarts('[N;a;r5;D2;h0;x1:1]:1:[N;r5;D2;h1;x2:3]:[N;r5;D2;h0;x1:2]:[C;r5]:[C;r5]:1')
    rules.append((q, 'tautomer_morgan_donor'))

    # === amidine: global H+bond resolution ===
    q = smarts('[N;h0,h1;z2:1]=[C;z2;x2,x3;!R:3]-[N;h1,h2;z1:2]')
    rules.append((q, 'amidine'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
