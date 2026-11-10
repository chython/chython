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
"""The three refusals of the polymer notations."""


__all__ = ['HelmError', 'BilnError', 'StereoMismatch']


class HelmError(ValueError):
    """A HELM string that does not parse, or a monomer graph HELM cannot spell."""


class BilnError(ValueError):
    """A BILN string that does not parse, or a molecule BILN cannot spell: inline SMILES, non-peptide chains."""


class StereoMismatch(ValueError):
    """A piece whose constitution is a library monomer and whose stereo is not."""
