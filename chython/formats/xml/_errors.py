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
"""The three ways an XML chemical record can end badly, kept apart because callers act on them
differently: :class:`MalformedXml` (not well-formed, or nothing to build a molecule from),
:class:`UnsupportedXml` (valid, states something this reader will not guess at) and
:class:`ForbiddenXml` (well-formed, refused by the entity/depth policy in :mod:`._tree`).
Everything short of these three is a log line."""


__all__ = ['XmlError', 'MalformedXml', 'UnsupportedXml', 'ForbiddenXml']


class XmlError(ValueError):
    """Base for all three, so a caller that does not care which can catch one thing."""


class MalformedXml(XmlError):
    """The document is not well-formed, or holds nothing a molecule can be built from."""


class UnsupportedXml(XmlError):
    """The document is valid and states a feature this reader refuses to guess at."""


class ForbiddenXml(XmlError):
    """The document is well-formed and the entity or depth policy declined to expand it.

    Distinct from :class:`MalformedXml` so a corpus scan can count refusals apart from damaged files,
    and so widening the policy (``allow_dtd=True``) does not widen what counts as a parse failure.
    """
