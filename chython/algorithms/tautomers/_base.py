# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ...periodictable import ListElement


def _stripped_rules():
    from ...containers import QueryContainer
    rules = []

    # Oxo-acid salts. [O,S,Se-]-[P,S,Cl,Se,Br,I]=O
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1)
    q.add_atom(ListElement(['P', 'S', 'Cl', 'Se', 'Br', 'I']))
    q.add_atom('O')
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 2)
    rules.append(q)

    # Phenole salts, alcoholates, carbon acids, nitrite. [O,S,Se-]-Ar
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1)
    q.add_atom(ListElement(['C', 'N']))
    q.add_bond(1, 2, 1)
    rules.append(q)

    # Nitrate. [O-]-[N+](=O)[O-]
    q = QueryContainer()
    q.add_atom('O', charge=-1)
    q.add_atom('N', charge=1)
    q.add_atom('O', charge=-1)
    q.add_atom('O')
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 1)
    q.add_bond(2, 4, 2)
    rules.append(q)

    # ions
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'F', 'Cl', 'Br', 'I']), charge=-1, neighbors=0)
    rules.append(q)
    return rules


def _rules():
    from ...containers import QueryContainer
    rules = _stripped_rules()

    # Guanidine
    q = QueryContainer()
    q.add_atom('N', heteroatoms=0)
    q.add_atom('C')
    q.add_atom('N')
    q.add_atom('N')
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    q.add_bond(2, 4, 1)
    rules.append(q)

    # Oxo-guanidine, Amino-guanidine
    q = QueryContainer()
    q.add_atom('N')
    q.add_atom('C')
    q.add_atom('N')
    q.add_atom('N')
    q.add_atom(ListElement(['O', 'N']))
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    q.add_bond(2, 4, 1)
    q.add_bond(1, 5, 1)
    rules.append(q)

    # O-alkyl-isourea, S-alkyl-isothiaurea
    q = QueryContainer()
    q.add_atom('N', heteroatoms=0)
    q.add_atom('C')
    q.add_atom('N')
    q.add_atom(ListElement(['O', 'S']), neighbors=2, hybridization=1, heteroatoms=0)
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    q.add_bond(2, 4, 1)
    rules.append(q)

    # Dialkyl imidocarbonate
    q = QueryContainer()
    q.add_atom('N', heteroatoms=0)
    q.add_atom('C')
    q.add_atom('O', neighbors=2, heteroatoms=0)
    q.add_atom('O', neighbors=2, heteroatoms=0)
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    q.add_bond(2, 4, 1)
    rules.append(q)

    # Amidine
    q = QueryContainer()
    q.add_atom('N', heteroatoms=0)
    q.add_atom('C', heteroatoms=2)
    q.add_atom('N')
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    rules.append(q)

    # O-alkyl-imidate (oxazoline)
    q = QueryContainer()
    q.add_atom('N', heteroatoms=0)
    q.add_atom('C', heteroatoms=2)
    q.add_atom('O', neighbors=2, heteroatoms=0)
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    rules.append(q)

    # Amidoxime. O-N=C([C,H])N
    q = QueryContainer()
    q.add_atom('N')
    q.add_atom('C', heteroatoms=2)
    q.add_atom('N')
    q.add_atom('O')
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    q.add_bond(1, 4, 1)
    rules.append(q)

    # Oxime, Hydrazone. [O,N]-N=C([C,H])[C,H]
    q = QueryContainer()
    q.add_atom('N')
    q.add_atom('C', heteroatoms=1, hybridization=2)
    q.add_atom(ListElement(['N', 'O']))
    q.add_bond(1, 2, 2)
    q.add_bond(1, 3, 1)
    rules.append(q)

    # Imine. [C,H]N=C([C,H])[C,H]
    q = QueryContainer()
    q.add_atom('N', heteroatoms=0)
    q.add_atom('C', heteroatoms=1, hybridization=2)
    q.add_bond(1, 2, 2)
    rules.append(q)

    # Alkyl amine, Hydroxylamine, Hydrazine
    q = QueryContainer()
    q.add_atom('N', neighbors=1)
    q.add_atom(ListElement(['C', 'O', 'N']), hybridization=1, heteroatoms=1)
    q.add_bond(1, 2, 1)
    rules.append(q)

    # Dialkyl amine, Alkyl hydroxylamine, Alkyl hydrazine
    q = QueryContainer()
    q.add_atom('N', neighbors=2)
    q.add_atom(ListElement(['C', 'O', 'N']), hybridization=1, heteroatoms=1)
    q.add_atom('C', hybridization=1, heteroatoms=1)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, 1)
    rules.append(q)

    # Trialkyl amine, Dialkyl-hydroxylamine, Dialkyl-hydrazine
    q = QueryContainer()
    q.add_atom('N')
    q.add_atom(ListElement(['C', 'O', 'N']), hybridization=1, heteroatoms=1)
    q.add_atom('C', hybridization=1, heteroatoms=1)
    q.add_atom('C', hybridization=1, heteroatoms=1)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, 1)
    q.add_bond(1, 4, 1)
    rules.append(q)

    # Pyridine. Imidazole. Triazole. :N:
    q = QueryContainer()
    q.add_atom('N', hybridization=4, hydrogens=0, neighbors=2)
    rules.append(q)
    return rules


stripped_rules = Proxy(_stripped_rules)
rules = Proxy(_rules)


__all__ = ['stripped_rules', 'rules']
