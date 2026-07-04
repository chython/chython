# -*- coding: utf-8 -*-
#
#  Copyright 2024-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
    from ... import smarts

    rules = {}

    # carbohydrides
    rules['terminal_alkene'] = smarts('[C;z2;x0;D1:1]=[C;z2;x0;D2,D3:2]')
    rules['alkene'] = smarts('[C;z2;x0;D2,D3:1]=[C;z2;x0;D2,D3:2]')
    rules['terminal_alkyne'] = smarts('[C;z3;x0;D1:1]#[C;x0;D2:2]')
    rules['alkyne'] = smarts('[C;z3;x0;D2:1]#[C;x0;D2:2]')

    # enamines/enol ethers: C=C bonded to N or O (for hydrogenation)
    rules['enamine'] = smarts('[C;z2;D2,D3:1]=[C;z2;D2,D3:2]-[N;z1;D2,D3:3]')
    rules['enol_ether'] = smarts('[C;z2;D2,D3:1]=[C;z2;D2,D3:2]-[O;D2:3]')

    rules['vicinal_diol'] = smarts('[O;D1;z1;x0:1]-[C;z1;x1]-[C;z1;x1]-[O;D1;z1;x0:2]')

    rules['benzylic_ch'] = smarts('[C;z1;h1,h2;D2,D3:1]-[C;a:2]')

    # halides
    rules['aryl_fluoride'] = smarts('[F;D1:100]-[C;a:1]')
    rules['aryl_chloride'] = smarts('[Cl;D1:100]-[C;a:1]')
    rules['aryl_bromide'] = smarts('[Br;D1:100]-[C;a:1]')
    rules['aryl_iodide'] = smarts('[I;D1:100]-[C;a:1]')

    rules['alkyl_fluoride'] = smarts('[F;D1:100][C;z1;x1:1]')
    rules['alkyl_chloride'] = smarts('[Cl;D1:100][C;z1;x1:1]')
    rules['alkyl_bromide'] = smarts('[Br;D1:100][C;z1;x1:1]')
    rules['alkyl_iodide'] = smarts('[I;D1:100][C;z1;x1:1]')

    rules['alkenyl_fluoride'] = smarts('[F;D1:100][C;z2;x1:1]=[C:2]')
    rules['alkenyl_chloride'] = smarts('[Cl;D1:100][C;z2;x1:1]=[C:2]')
    rules['alkenyl_bromide'] = smarts('[Br;D1:100][C;z2;x1:1]=[C:2]')
    rules['alkenyl_iodide'] = smarts('[I;D1:100][C;z2;x1:1]=[C:2]')

    rules['alkynyl_fluoride'] = smarts('[F;D1:100][C;z3;x1:1]')
    rules['alkynyl_chloride'] = smarts('[Cl;D1:100][C;z3;x1:1]')
    rules['alkynyl_bromide'] = smarts('[Br;D1:100][C;z3;x1:1]')
    rules['alkynyl_iodide'] = smarts('[I;D1:100][C;z3;x1:1]')

    # pseudohalides
    rules['aryl_triflate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;a:1])-[C;D4](F)(F)F')
    rules['aryl_mesylate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;a:1])-[C;D1]')
    rules['aryl_tosylate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;a:1])-[C;a]:1:[C;D2]:[C;D2]:[C](-[C;D1]):[C;D2]:[C;D2]:1')

    rules['alkyl_triflate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;z1;x1:1])-[C;D4](F)(F)F')
    rules['alkyl_mesylate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;z1;x1:1])-[C;D1]')
    rules['alkyl_tosylate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;z1;x1:1])-[C;a]:1:[C;D2]:[C;D2]:[C](-[C;D1]):[C;D2]:[C;D2]:1')

    # boronic acids and esters
    rules['aryl_boronic_acid'] = smarts('[B;D3;z1;x2:100](-[O;D1])(-[O;D1])-;!@[C;a:1]')
    rules['aryl_boronic_ester'] = smarts('[B;D3;z1;x2:100](-[O;D2;x1])(-[O;D2;x1])-;!@[C;a:1]')

    rules['alkyl_boronic_acid'] = smarts('[B;D3;z1;x2:100](-[O;D1])(-[O;D1])-;!@[C;z1;x1:1]')
    rules['alkyl_boronic_ester'] = smarts('[B;D3;z1;x2:100](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z1;x1:1]')

    rules['alkenyl_boronic_acid'] = smarts('[B;D3;z1;x2:100](-[O;D1])(-[O;D1])-;!@[C;z2;x1:1]=[C:2]')
    rules['alkenyl_boronic_ester'] = smarts('[B;D3;z1;x2:100](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z2;x1:1]=[C:2]')

    rules['alkynyl_boronic_acid'] = smarts('[B;D3;z1;x2:100](-[O;D1])(-[O;D1])-;!@[C;z3;x1:1]')
    rules['alkynyl_boronic_ester'] = smarts('[B;D3;z1;x2:100](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z3;x1:1]')

    # molander salts (trifluoroborates)
    rules['aryl_molander_salt'] = smarts('[B;D4;z1;x3;-:100](F)(F)(F)-;!@[C;a:1]')
    rules['alkyl_molander_salt'] = smarts('[B;D4;z1;x3;-:100](F)(F)(F)-;!@[C;z1:1]')
    rules['alkenyl_molander_salt'] = smarts('[B;D4;z1;x3;-:100](F)(F)(F)-;!@[C;z2:1]=[C:2]')
    rules['alkynyl_molander_salt'] = smarts('[B;D4;z1;x3;-:100](F)(F)(F)-;!@[C;z3;x1:1]')

    # alcohols and phenols
    rules['primary_alcohol'] = smarts('[O;D1;z1;x0:1][C;D2;x1;z1:2]')
    rules['secondary_alcohol'] = smarts('[O;D1;z1;x0:1][C;D3;x1;z1:2]')
    rules['tertiary_alcohol'] = smarts('[O;D1;z1;x0:1][C;D4;x1;z1:2]')
    # tertiary alcohol with adjacent sp3 CH (eliminable): for dehydration
    rules['tertiary_alcohol_with_alpha_h'] = smarts('[O;D1;z1;x0:1]-[C;D4;x1;z1:2]-[C;z1;h1,h2:3]')
    rules['phenol'] = smarts('[O;D1;z1;x0:1]-[C;a:2]')

    # aldehydes and ketones
    rules['aldehyde'] = smarts('[O;z2;x0:2]=[C;D2;x1;z2:1]')
    rules['ketone'] = smarts('[O;z2;x0:2]=[C;D3;x1;z2:1]')
    # enal: alpha,beta-unsaturated aldehyde (for Doebner-Miller)
    rules['enal'] = smarts('[O;z2;x0:2]=[C;D2;x1;z2:1]-[C;z2;x0:3]=[C;x0:4]')
    # fisher, friedlander
    rules['alpha_ketone'] = smarts('[O;z2;x0:2]=[C;D3;x1:1]-[C;z1;D1,D2;x0:3]')
    # hantzsch thiazole, imidazo[1,2-a]pyridine
    rules['alpha_haloketone'] = smarts('[O;z2;x0:2]=[C;D3;x1:1]-[C;z1;D2,D3;x1:3]([Cl,Br,I;D1:100])')
    # imidazo[1,2-a]pyridine from alpha-haloester (ester O masked, carbonyl O deleted)
    rules['alpha_haloester'] = smarts('[O;D2;x0;M]-[C;D3;x2;z2:1](=[O:2])-[C;z1;D2,D3;x1:3]([Cl,Br,I;D1:100])')

    rules['1_2_diketone'] = smarts('[O;z2;x0:2]=[C;D3;x1:1]-[C;z2;x1;D3:3]=[O:4]')
    rules['1_3_diketone'] = smarts('[O;z2;x0:2]=[C;D3;x1:1]-[C;z1;D2,D3:5]-[C;z2;x1;D3:3]=[O:4]')
    rules['1_4_diketone'] = smarts('[O;z2;x0:2]=[C;D3;x1:1]-[C;z1;D2,D3:5]-[C;z1;D2,D3:6]-[C;z2;x1;D3:3]=[O:4]')
    rules['beta_ketoester'] = smarts('[O;z2;x0:2]=[C;D3;x1:1]-[C;z1;D2;x0:5]-[C;z2;x2;D3:3](=[O:4])[O;D2:100]')

    # acids
    rules['alkyl_carboxylic_acid'] = smarts('[O;D1;z1;x0:100][C;z2;x2;D3:1](=[O:2])[C;z1:3]')
    rules['aryl_carboxylic_acid'] = smarts('[O;D1;z1;x0:100][C;z2;x2;D3:1](=[O:2])[C;a:3]')
    rules['carboxylic_acid'] = smarts('[O;D1;z1;x0:100][C;z2;x2;D3:1]=[O:2]')
    rules['acyl_chloride'] = smarts('[Cl:100][C;z2;x2;D3:1]=[O:2]')
    rules['acyl_fluoride'] = smarts('[F:100][C;z2;x2;D3:1]=[O:2]')
    rules['chloroformate'] = smarts('[Cl:100][C;z2;x3;D3:1](=[O:2])-[O;D2:3]')
    rules['fluoroformate'] = smarts('[F:100][C;z2;x3;D3:1](=[O:2])-[O;D2:3]')
    rules['carbamoyl_chloride'] = smarts('[Cl:100][C;z2;x3;D3:1](=[O:2])-[N;D2,D3:3]')
    rules['carbamoyl_fluoride'] = smarts('[F:100][C;z2;x3;D3:1](=[O:2])-[N;D2,D3:3]')

    # amines
    rules['primary_amine'] = smarts('[N;D1;z1;x0:1][C;z1:2]')
    rules['primary_aniline'] = smarts('[N;D1;z1;x0:1][C;a:2]')
    # NH2 on sp2 C=N (e.g. pyrazoline, amidine-like)
    rules['primary_amidine_amine'] = smarts('[N;D1;z1;x0:1]-[C;z2:2]=[N;M]')
    rules['secondary_amine'] = smarts('[N;D2;z1;x0:1]([C;z1:2])[C;z1:3]')
    rules['secondary_aniline'] = smarts('[N;D2;z1;x0:1]([C;a:2])[C;z1:3]')
    rules['biaryl_aniline'] = smarts('[N;D2;z1;x0:1]([C;a:2])[C;a:3]')

    # esters and amides
    rules['ester'] = smarts('[O;z2;x0:2]=[C;D3;x2;z2:1]-[O;D2;x0:100]')
    rules['primary_amide'] = smarts('[N;D1;z1;x0:1][C;z2;x2;D3:2]=[O:3]')
    rules['secondary_amide'] = smarts('[N;D2;z1;x0:1][C;z2;x2;D3:2]=[O:3]')

    # sulfonyl
    rules['sulfonyl_chloride'] = smarts('[S;D4:1](=[O:2])(=[O:3])[Cl;D1:100]')
    rules['sulfonyl_fluoride'] = smarts('[S;D4:1](=[O:2])(=[O:3])[F;D1:100]')
    rules['sulfonamide'] = smarts('[S;x3;D4:1](=[O:2])(=[O:3])-[N;z1:100]')
    rules['sulfonyl_anhydride'] = smarts('[S;x3;D4:1](=[O:2])(=[O:3])-[O:100]-[S;x3;D4](=O)(=O)')

    # nitrogen functional groups
    rules['nitrile'] = smarts('[N;D1;z3;x0:2]#[C;D2;x1:1]')
    rules['azide'] = smarts('[N;x1;D2:1]=[N+:2]=[N-:3]')
    rules['isocyanate'] = smarts('[N;z2;x0;D2:1]=[C:2]=[O:3]')
    rules['isocyano'] = smarts('[N;D2;x0;+:1]#[C;-:2]')
    rules['guanidine'] = smarts('[N;z1;x0:1][C;!R:2]([N;z1;x0:3])=[N;x0:4]')
    rules['nitro'] = smarts('[N;D3;x2;+:1]([O;-:2])=[O:3]')

    # grignard reagents (RMgX)
    rules['alkyl_grignard'] = smarts('[Mg;D2:100](-[F,Cl,Br,I])-[C;z1:1]')
    rules['aryl_grignard'] = smarts('[Mg;D2:100](-[F,Cl,Br,I])-[C;a:1]')
    rules['alkenyl_grignard'] = smarts('[Mg;D2:100](-[F,Cl,Br,I])-[C;z2:1]=[C:2]')

    # organozinc reagents (RZnX)
    rules['alkyl_zinc'] = smarts('[Zn;D2:100](-[F,Cl,Br,I])-[C;z1:1]')
    rules['aryl_zinc'] = smarts('[Zn;D2:100](-[F,Cl,Br,I])-[C;a:1]')
    rules['alkenyl_zinc'] = smarts('[Zn;D2:100](-[F,Cl,Br,I])-[C;z2:1]=[C:2]')

    # boronate alkyl halides: B-CH2-X (one-carbon, for SN2)
    rules['boronate_alkyl_chloride'] = smarts('[Cl;D1:100]-[C;z1;D2;x2:1]-[B;M]')
    rules['boronate_alkyl_bromide'] = smarts('[Br;D1:100]-[C;z1;D2;x2:1]-[B;M]')
    rules['boronate_alkyl_iodide'] = smarts('[I;D1:100]-[C;z1;D2;x2:1]-[B;M]')

    # stannanes (R-SnR3)
    rules['aryl_stannane'] = smarts('[Sn;D4;z1:100]-;!@[C;a:1]')
    rules['alkenyl_stannane'] = smarts('[Sn;D4;z1:100]-;!@[C;z2:1]=[C:2]')
    rules['alkyl_stannane'] = smarts('[Sn;D4;z1:100]-;!@[C;z1;D2,D3,D4:1]')

    # silanes (R-SiR3, for Hiyama coupling)
    rules['aryl_silane'] = smarts('[Si;D4:100]-;!@[C;a:1]')
    rules['alkenyl_silane'] = smarts('[Si;D4:100]-;!@[C;z2:1]=[C:2]')
    rules['alkynyl_silane'] = smarts('[Si;D4:100]-;!@[C;z3:1]#[C:2]')

    # phosphorus ylides and phosphonates
    rules['phosphonium_ylide'] = smarts('[P;D4;z2;x0:100]=[C:1]')
    rules['phosphonate'] = smarts('[P;D4;x3:100](=O)([O;D2;x1])([O;D2;x1])-[C:1]')

    # weinreb amide
    rules['weinreb_amide'] = smarts('[O:2]=[C;D3;x2:1]-[N;D3;x1:100][O;D2;x1]')

    # arene C-H (for electrophilic aromatic substitution)
    rules['arene_ch'] = smarts('[C;a;D2:1]')

    # sulfur
    rules['thiol'] = smarts('[S;x0;D1;z1:1][C;z1:2]')
    rules['thioether'] = smarts('[S;D2;z1;x0:1]([C:2])[C:3]')
    rules['sulfoxide'] = smarts('[S;D3;z2:1](=[O:2])([C:3])[C:4]')
    rules['sulfone'] = smarts('[S;D4:1](=[O:2])(=[O:3])([C:4])[C:5]')

    # lactam halides: ring-closure patterns for activated C(sp2)-X in 6-membered lactam rings
    # covers pyridone, pyridazinone, quinazolinone, pyrimidinone, benzo-fused lactams
    # 6 positional patterns per halide (all relative positions of X to N-C=O in ring)
    for _x, _name in (('F', 'fluoride'), ('Cl', 'chloride'), ('Br', 'bromide'), ('I', 'iodide')):
        # X=C-C/N-C(=O)-N-C/N  e.g. CN1C=CC(Br)=CC1=O
        rules[f'lactam_1_{_name}'] = smarts(f'[{_x};D1:100]-[C;z2;r6:1]1=[C,N;z2;M]-[C;D3;z2;M]-[N;D3;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]1')
        # X=C-N-C(=O)-C/N-C/N  e.g. CN1C(Br)=CC=CC1=O
        rules[f'lactam_2_{_name}'] = smarts(f'[{_x};D1:100]-[C;z2;r6:1]=1-[N;D3;M]-[C;D3;z2;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]-[C,N;z2;M]=1')
        # X=C-N-C/N-C/N-C(=O)  e.g. CN1C=CC(=O)C=C1Br
        rules[f'lactam_3_{_name}'] = smarts(f'[{_x};D1:100]-[C;z2;r6:1]=1-[N;D3;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]-[C;D3;z2;M]-[C,N;z2;M]=1')
        # X=C-C(=O)-N-C/N-C/N  e.g. CN1C=CC=C(Br)C1=O
        rules[f'lactam_4_{_name}'] = smarts(f'[{_x};D1:100]-[C;z2;r6:1]=1-[C;D3;z2;M]-[N;D3;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]-[C,N;z2;M]=1')
        # X=C-C/N-C/N-N-C(=O)  e.g. CN1C=CC(=O)C(Br)=C1
        rules[f'lactam_5_{_name}'] = smarts(f'[{_x};D1:100]-[C;z2;r6:1]1=[C,N;z2;M]-[N;D3;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]-[C;D3;z2;M]1')
        # X=C-C/N-N-C(=O)-C/N  e.g. CN1C=C(Br)C=CC1=O
        rules[f'lactam_6_{_name}'] = smarts(f'[{_x};D1:100]-[C;z2;r6:1]1=[C,N;z2;M]-[N;D3;M]-[C;D3;z2;M]-[C,N;z2,z4;M]=,:[C,N;z2,z4;M]1')

    # pyrrole. for tautomerism handling H not in template.
    rules['pyrrole'] = smarts('[N;h1;D2;a;r5:1]')
    rules['pyrazole'] = smarts('[N;h1;D2;a;r5:1]:[N;h0;D2;r5:2]')
    rules['imidazole'] = smarts('[N;h1;D2;a;r5:1]:[A:2]:[N;h0;D2;r5:3]')

    # heterocycles (for screening)
    rules['isoxazole'] = smarts('[O;a;D2;r5:1]:[N;a;D2;r5:2]')
    rules['pyridazine'] = smarts('[N;a;D2;r6:1]:[N;a;D2;r6:2]')

    # hydrazines
    rules['alkyl_hydrazine'] = smarts('[N;D1;z1;x1:2]-[N;D2;z1;x1:1]-[C;z1:3]')
    rules['aryl_hydrazine'] = smarts('[N;D1;z1;x1:2]-[N;D2;z1;x1:1]-[C;a:3]:[C;a;D2:4]')

    # hydrazone (C=N-NH-R, product of carbonyl + hydrazine condensation)
    rules['hydrazone'] = smarts('[C;z2:1]=[N;D2;z2;x1:2]-[N;D2;z1;x1:3]')

    # sulfonylhydrazone (C=N-NH-SO2R, e.g. tosylhydrazone; precursor to diazo via Bamford-Stevens)
    # N:3 bonded to N and S -> x2
    rules['sulfonylhydrazone'] = smarts('[C;z2:1]=[N;D2;z2;x1:2]-[N;D2;z1;x2:3]-[S;D4;x3:100](=[O])(=[O])')

    # thioamide (for Hantzsch thiazole)
    rules['thioamide'] = smarts('[S;z2;x0;D1:2]=[C;D3;x2:1]-[N;D1:3]')

    # ortho-bifunctional arenes
    rules['o_diaminoarene'] = smarts('[N;D1;z1;x0:1]-[C;a:3]:[C;a:4]-[N;D1,D2;z1;x0:2]')
    rules['o_aminophenol'] = smarts('[N;D1;z1;x0:1]-[C;a:3]:[C;a:4]-[O;D1:2]')
    rules['o_aminothiophenol'] = smarts('[N;D1;z1;x0:1]-[C;a:3]:[C;a:4]-[S;D1:2]')
    rules['o_aminobenzaldehyde'] = smarts('[N;D1;z1;x0:1]-[C;a:3]:[C;a:4]-[C;D2;z2;x1:5]=[O:6]')
    rules['anthranilic_acid'] = smarts('[N;D1;z1;x0:1]-[C;a:3]:[C;a:4]-[C;z2;x2;D3:5](=[O:6])[O;D1:100]')

    # amidoxime: RC(=NH)NHOH canonical form (for 1,2,4-oxadiazole)
    rules['amidoxime'] = smarts('[N;D1;z2;x0:3]=[C;D3;x2:1]-[N;D2;z1;x1:2]-[O;D1:4]')

    # amidine: RC(=NH)NH2 (for pyrimidine)
    rules['amidine'] = smarts('[N;D1;z1;x0:3]-[C;D3;z2;x2:1]=[N;D1:2]')
    # urea/thiourea (for Biginelli)
    rules['urea'] = smarts('[N;D1;z1;x0:1]-[C;D3;z2;x3:2](=[O:3])-[N;D1:4]')
    rules['thiourea'] = smarts('[N;D1;z1;x0:1]-[C;D3;z2;x3:2](=[S;D1:3])-[N;D1:4]')

    # beta-arylethylamine (for Pictet-Spengler)
    rules['beta_arylethylamine'] = smarts('[N;D1;z1;x0:1]-[C;z1:2]-[C;z1:3]-[C;a:4]:[C;a;D2:5]')
    # 2-aminopyridine / 2-aminoazine (for GBB, imidazo[1,2-a]pyridine)
    rules['aminopyridine'] = smarts('[N;D1;z1;x0:1]-[C;a:2]:[N;a;h0;D2:3]')

    # amino alcohol (for oxazoline formation): H2N-C-C-OH
    rules['amino_alcohol'] = smarts('[N;D1;z1;x0:1]-[C;z1:2]-[C;z1:3]-[O;D1:4]')

    # hydroxamic acid: R-C(=O)-NH-OH
    rules['hydroxamic_acid'] = smarts('[O;D1;z1;x1:1]-[N;D2;z1;x1:2]-[C;z2;x2:3]=[O:4]')

    # oxime: C=N-OH
    rules['oxime'] = smarts('[O;D1;z1;x1:1]-[N;D2;z2;x1:2]=[C:3]')

    # O-alkylhydroxylamine: R-O-NH2 (for oxime ether formation)
    rules['O_alkylhydroxylamine'] = smarts('[N;D1;z1;x1:1]-[O;D2;z1;x1:2]-[C:3]')

    # alpha-isocyano (for Van Leusen oxazole)
    rules['tosyl_isocyanide'] = smarts('[C;-;D1:2]#[N;+;D2:1]-[C;D2,D3;z1;x2:3]-[S;D4;x2:100](=O)=O')

    # thioester (for Liebeskind-Srogl): R-C(=O)-S-R'
    rules['thioester'] = smarts('[O;z2;x0:2]=[C;D3;x2;z2:1]-[S;D2;z1;x0:100]')

    # active methylene (for Knoevenagel): CH flanked by 2 EWGs
    rules['active_methylene'] = smarts('[C;z1;D2,D3;x0:1](-[C;z2,z3;x1,x2:2])-[C;z2,z3;x1,x2:3]')

    # aniline with ortho C-H (for Doebner-Miller)
    rules['aniline_ortho_ch'] = smarts('[N;D1;z1;x0:1]-[C;a:2]:[C;a;D2:3]')

    # ortho-haloaniline (for Larock indole)
    rules['o_haloaniline'] = smarts('[N;D1;z1;x0:1]-[C;a:2]:[C;a:3]-[Cl,Br,I;D1:100]')

    # tertiary amine (for N-oxidation): NR3 where all substituents are C
    rules['tertiary_amine'] = smarts('[N;D3;z1;x0:1]')

    # pyridine-like nitrogen (for N-oxidation): aromatic N with no H, degree 2
    rules['pyridine_n'] = smarts('[N;a;D2;h0:1]')

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
