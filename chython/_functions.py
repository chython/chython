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


__all__ = ['lazy_product']
