# -*- coding: utf-8 -*-
#
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
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
from .acceptor import queries as acceptor_queries, banned as acceptor_banned
from .acidic import queries as acidic_queries
from .aromatic import queries as aromatic_queries
from .basic import queries as basic_queries, banned as basic_banned
from .donor import queries as donor_queries
from .halogen import queries as halogen_queries


__all__ = [
    'acceptor_queries', 'acceptor_banned',
    'acidic_queries',
    'aromatic_queries',
    'basic_queries', 'basic_banned',
    'donor_queries',
    'halogen_queries'
]
