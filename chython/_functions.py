# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import product
from sys import version_info


# lazy itertools.product with diagonal combination precedence
def lazy_product(*args):
    if len(args) == 1:
        for x in args[0]:
            yield x,
    elif not args:
        yield ()
    else:
        gens = [iter(x) for x in args]
        empty = [False] * len(args)
        pools = [[] for _ in range(len(args))]
        indices = set()

        reached = 0
        while True:
            out = []
            ind = []
            for n, (p, g, e) in enumerate(zip(pools, gens, empty)):
                if e:
                    out.append(p[-1])
                else:
                    try:
                        x = next(g)
                    except StopIteration:
                        if not p:  # one of gens empty
                            return
                        reached += 1
                        if reached == len(args):
                            break
                        out.append(p[-1])
                        empty[n] = True
                    else:
                        p.append(x)
                        out.append(x)
                ind.append(len(p) - 1)
            else:
                yield tuple(out)
                indices.add(tuple(ind))
                continue
            break

        for ind in product(*(range(len(p)) for p in pools)):
            if ind in indices:
                continue
            yield tuple(p[x] for x, p in zip(ind, pools))


if version_info[1] >= 8:
    tuple_hash = hash
else:
    def tuple_hash(v) -> int:
        """
        Python 3.8 hash for tuples implemented on python.
        """
        acc = 0x27D4EB2F165667C5
        for el in v:
            if isinstance(el, tuple):
                lane = tuple_hash(el)
            else:
                lane = hash(el)
            if lane == -1:
                return -1
            elif lane < 0:
                lane += 0x10000000000000000  # to unsigned

            acc += lane * 0xC2B2AE3D27D4EB4F
            acc %= 0x10000000000000000
            acc = (acc << 31) % 0x10000000000000000 | (acc >> 33)
            acc *= 0x9E3779B185EBCA87

        acc += len(v) ^ 2870177450013471926  # 0xC2B2AE3D27D4EB4F ^ 3527539
        acc %= 0x10000000000000000

        if acc == 0xFFFFFFFFFFFFFFFF:
            return 1546275796
        elif acc > 0x7FFFFFFFFFFFFFFF:  # (1 << 63) - 1 = largest positive number
            return acc - 0x10000000000000000
        return acc


__all__ = ['lazy_product', 'tuple_hash']
