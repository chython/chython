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


def _make_reactor(rxn_name, fg_name, output_fg, product_smarts):
    from ._functional import rules
    from ...reactor import Reactor
    from ... import smarts

    q = rules[fg_name]
    product = smarts(product_smarts)
    return rxn_name, fg_name, output_fg, Reactor((q,), (product,),
                                                  delete_atoms=True, one_shot=True, fix_broken_pyrroles=True)


def _rules():
    rules = []

    # primary_alcohol → aldehyde (Swern, Dess-Martin, PCC)
    rules.append(_make_reactor('alcohol_to_aldehyde', 'primary_alcohol', 'aldehyde', '[A:1]=[A:2]'))

    # secondary_alcohol → ketone (Dess-Martin, Jones, PCC)
    rules.append(_make_reactor('alcohol_to_ketone', 'secondary_alcohol', 'ketone', '[A:1]=[A:2]'))

    # aldehyde → carboxylic_acid (Pinnick, Jones, KMnO4)
    rules.append(_make_reactor('aldehyde_to_acid', 'aldehyde', 'carboxylic_acid', '[A:2]=[A:1]-[O:20]'))

    # alkene → 1,2-diol (Sharpless dihydroxylation, OsO4/KMnO4)
    rules.append(_make_reactor('dihydroxylation', 'terminal_alkene', 'vicinal_diol', '[A:1](-[O:20])-[A:2](-[O:21])'))
    rules.append(_make_reactor('dihydroxylation', 'alkene', 'vicinal_diol', '[A:1](-[O:20])-[A:2](-[O:21])'))

    # thioether → sulfoxide (mCPBA, H2O2, NaIO4)
    rules.append(_make_reactor('thioether_to_sulfoxide', 'thioether', 'sulfoxide', '[A:1](=[O:20])(-[A:2])-[A:3]'))

    # thioether → sulfone (excess mCPBA, H2O2/AcOH, Oxone)
    rules.append(_make_reactor('thioether_to_sulfone', 'thioether', 'sulfone', '[A:1](=[O:20])(=[O:21])(-[A:2])-[A:3]'))

    # sulfoxide → sulfone
    rules.append(_make_reactor('sulfoxide_to_sulfone', 'sulfoxide', 'sulfone', '[A:1](=[A:2])(=[O:20])(-[A:3])-[A:4]'))

    # tertiary amine → N-oxide (mCPBA, H2O2)
    rules.append(_make_reactor('nitrogen_oxidation', 'tertiary_amine', None, '[N;+:1]([O;-:20])'))

    # pyridine → N-oxide (mCPBA, H2O2)
    rules.append(_make_reactor('nitrogen_oxidation', 'pyridine_n', None, '[N;+:1]([O;-:20])'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
