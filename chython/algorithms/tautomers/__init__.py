# -*- coding: utf-8 -*-
#
#  Copyright 2020-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Nail Samikaev <samikaevn@yandex.ru>
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
from collections import deque
from typing import TYPE_CHECKING, Iterator, Union
from .acid_base import *
from .heteroarenes import *
from .keto_enol import *


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Tautomers(AcidBase, HeteroArenes, KetoEnol):
    """
    Oxides and sulphides ignored.
    """
    __slots__ = ()

    def enumerate_tautomers(self: Union['MoleculeContainer', 'Tautomers'], *, prepare_molecules=True, zwitter=True,
                            partial=False, increase_aromaticity=True, keep_sugars=True, heteroarenes=True,
                            keto_enol=True, limit: int = 1000) -> Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule.

        :param prepare_molecules: Standardize structures for correct processing
        :param zwitter: Do zwitter-ions enumeration
        :param partial: Allow OC=CC=C>>O=CCC=C or O=CC=CC>>OC=C=CC
        :param increase_aromaticity: prevent aromatic ring destruction
        :param keep_sugars: prevent carbonyl moving in sugars
        :param heteroarenes: enumerate heteroarenes
        :param keto_enol: enumerate keto-enols
        :param limit: Maximum attempts count
        """
        if limit < 1:
            raise ValueError('limit should be greater or equal 1')
        counter = 0

        copy = self.copy(keep_sssr=True, keep_components=True)
        if prepare_molecules:  # transform to kekule form without hydrogens
            copy.kekule()
            copy.implicify_hydrogens(_fix_stereo=False)

        # transform to thiele to prevent duplicates and dearomatization
        thiele = copy.copy(keep_sssr=True, keep_components=True)
        thiele.thiele(fix_tautomers=False)
        yield thiele  # return original structure as first tautomer

        seen = {thiele: None}  # value is parent molecule - required for preventing migrations in sugars.

        # first try to neutralize
        if copy.neutralize(_fix_stereo=False):  # found neutral form
            thiele = copy.copy(keep_sssr=True, keep_components=True)
            thiele.thiele(fix_tautomers=False)
            yield thiele
            counter += 1
            seen[thiele] = None

        # lets iteratively do keto-enol transformations.
        rings_count = len(thiele.aromatic_rings)  # increase rings strategy.
        if keto_enol:
            queue = deque([(copy, thiele)])
        else:
            queue = None
        new_queue = [thiele]  # new_queue - molecules suitable for hetero-arenes enumeration.
        # store aromatic form to seen. kekule forms not suitable for duplicate checking.

        while queue:
            current, thiele_current = queue.popleft()
            for mol, ket in current._enumerate_keto_enol_tautomers(partial):
                thiele = mol.copy(keep_sssr=True, keep_components=True)
                thiele.thiele(fix_tautomers=False)

                if thiele not in seen:
                    seen[thiele] = current
                    rc = len(thiele.aromatic_rings)
                    if increase_aromaticity:
                        if rc < rings_count:  # skip aromatic rings destruction
                            continue
                        elif rc > rings_count:  # higher aromaticity found. flush old queues.
                            rings_count = rc
                            queue = deque([(mol, thiele)])
                            new_queue = [thiele]
                            copy = mol  # new entry point.
                            yield thiele
                            break
                    if keep_sugars and current is not copy and ket:
                        # prevent carbonyl migration in sugars. skip entry point.
                        # search alpha hydroxy ketone inversion
                        before = seen[thiele_current]._sugar_groups
                        if any((k, e) in before for e, k in mol._sugar_groups):
                            continue

                    queue.append((mol, thiele))
                    new_queue.append(thiele)
                    yield thiele
                counter += 1
                if counter == limit:
                    return

        # heteroarenes tautomers enumeration
        if heteroarenes:
            queue = deque(new_queue)
        while queue:
            current = queue.popleft()
            for mol in current._enumerate_hetero_arene_tautomers():
                if mol not in seen:
                    seen[mol] = None
                    queue.append(mol)
                    new_queue.append(mol)  # new hetero-arenes also should be included to this list.
                    yield mol
                counter += 1
                if counter == limit:
                    return

        # zwitter-ions enumeration
        if zwitter:
            queue = deque(new_queue)
        while queue:
            current = queue.popleft()
            for mol in current._enumerate_zwitter_tautomers():
                if mol not in seen:
                    seen[mol] = None
                    queue.append(mol)
                    yield mol
                counter += 1
                if counter == limit:
                    return

    def enumerate_charged_tautomers(self: 'MoleculeContainer', *, prepare_molecules=True, partial=False,
                                    increase_aromaticity=True, keep_sugars=True, heteroarenes=True,
                                    keto_enol=True, deep: int = 4, limit: int = 1000):
        """
        Enumerate tautomers and protonated-deprotonated forms.
        Better to use on neutralized non-ionic molecules.

        See `enumerate_tautomers` and `enumerate_charged_forms` params description.
        """
        count = 0
        for t in self.enumerate_tautomers(prepare_molecules=prepare_molecules, zwitter=False, partial=partial,
                                          increase_aromaticity=increase_aromaticity, keep_sugars=keep_sugars,
                                          heteroarenes=heteroarenes, keto_enol=keto_enol, limit=limit):
            yield t
            count += 1
            if count == limit:
                return
            for c in t.enumerate_charged_forms(deep=deep, limit=limit):
                yield c
                count += 1
                if count == limit:
                    return


__all__ = ['Tautomers']
