# -*- coding: utf-8 -*-
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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


def _stripped_rules():
    from ... import smarts

    rules = []

    # Oxo-acid salts
    q = smarts('[O,S,Se;D1;z1;-][C,Si,N,P,S,Se,Cl,Br,I]=O')
    rules.append(q)

    # Thiophosphate salts
    q = smarts('P([O,S,Se;D1;-:1])=[S,Se]')
    rules.append(q)

    # Phenole salts, alcoholates
    q = smarts('[O,S,Se;D1;z1;-][C,N]')
    rules.append(q)

    # Nitrate
    q = smarts('[N;D3;z2;+]([O;-:1])([O-])=O')
    rules.append(q)

    # imide
    q = smarts('[N;D1,D2;z1;-]')
    rules.append(q)

    # ions
    q = smarts('[O,S,F,Cl,Br,I;D0;-]')
    rules.append(q)
    return rules


def _rules():
    from ... import smarts

    rules = _stripped_rules()

    # Guanidine
    q = smarts('[N;x0;z2]=C(N)N')
    rules.append(q)

    # Oxo-guanidine, Amino-guanidine
    q = smarts('[N;D2;x1;z2]([O,N])=C(N)N')
    rules.append(q)

    # O-alkyl-isourea, S-alkyl-isothiaurea
    q = smarts('[N;x0;z2]=C([O,S;D2;x0;z1])N')
    rules.append(q)

    # Dialkyl imidocarbonate
    q = smarts('[N;x0;z2]=C([O;D2;x0])[O;D2;x0]')
    rules.append(q)

    # Amidine
    q = smarts('[N;x0;z2]=[C;x2]N')
    rules.append(q)

    # O-alkyl-imidate (oxazoline)
    q = smarts('[N;x0;z2]=[C;x2][O;D2;x0]')
    rules.append(q)

    # Amidoxime. O-N=C([C,H])N
    q = smarts('[N;D2;x1;z2](O)=[C;x2]N')
    rules.append(q)

    # Oxime, Hydrazone. [O,N]-N=C([C,H])[C,H]
    q = smarts('[N;D2;x1;z2]([N,O])=[C;x1;z2]')
    rules.append(q)

    # Imine
    q = smarts('[N;x0;z2]=[C;x1;z2]')
    rules.append(q)

    # Alkyl amine, Hydroxylamine, Hydrazine
    q = smarts('[N;D1;z1][C,N,O;x1;z1]')
    rules.append(q)

    # Dialkyl amine, Alkyl hydroxylamine, Alkyl hydrazine
    q = smarts('[N;D2;z1]([C,N,O;x1;z1])[C;x1;z1]')
    rules.append(q)

    # Trialkyl amine, Dialkyl-hydroxylamine, Dialkyl-hydrazine
    q = smarts('[N;D3;z1]([C,N,O;x1;z1])([C;x1;z1])[C;x1;z1]')
    rules.append(q)

    # Pyridine. Imidazole. Triazole. :N:
    q = smarts('[N;a;h0;D2]')
    rules.append(q)
    return rules


stripped_rules = Proxy(_stripped_rules)
rules = Proxy(_rules)


__all__ = ['stripped_rules', 'rules']
