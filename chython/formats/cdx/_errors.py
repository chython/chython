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
"""The two ways a CDX or CDXML document is refused, mirroring the CTfile pair."""


__all__ = ['CdxError', 'MalformedCdx', 'UnsupportedCdx']


class CdxError(ValueError):
    """Base for both, so a caller that does not care which can catch one thing."""


class MalformedCdx(CdxError):
    """The document is damaged beyond what the recovery rules cover: no header, no structure at all."""


class UnsupportedCdx(CdxError):
    """The document is valid and states a query feature: a structure container cannot hold a constraint."""
