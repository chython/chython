# -*- coding: utf-8 -*-
#
#  Copyright 2019-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""Vendored from ``CachedMethods`` to drop a third-party dep that broke
chython on PyPI release boundaries (priority-rules bug #5).

No locks — dict get/setitem are atomic on both GIL and free-threaded
3.14+, wrapped functions are pure, duplicate computation on cache miss
is benign. Lock-free keeps cached state picklable.
"""
from collections.abc import Mapping
from functools import wraps


_SENTINEL = object()


def _freeze(value):
    """Freeze mutable returns so cached values can't be silently mutated."""
    if isinstance(value, list):
        return tuple(value)
    elif isinstance(value, set):
        return frozenset(value)
    elif isinstance(value, dict):
        return FrozenDict(value)
    return value


class FrozenDict(Mapping):
    """Immutable dict — returned by ``_freeze``."""
    __slots__ = ('__d',)

    def __init__(self, *args, **kwargs):
        self.__d = dict(*args, **kwargs)

    def __iter__(self):
        return iter(self.__d)

    def __len__(self):
        return len(self.__d)

    def __getitem__(self, key):
        return self.__d[key]

    def __repr__(self):
        return repr(self.__d)

    def copy(self):
        """Return a mutable shallow copy as a regular ``dict``."""
        return self.__d.copy()


class cached_property:
    """Per-instance cached property. Requires ``__dict__`` on the host.
    For slotted classes, use ``class_cached_property``."""

    def __init__(self, func):
        self.__doc__ = getattr(func, "__doc__")
        self.func = func
        name = func.__name__
        if name.startswith('__') and not name.endswith('__'):
            name = f'_{func.__qualname__.split(".")[-2]}{name}'
        self.name = name

    def __get__(self, obj, cls):
        if obj is None:
            return self
        value = obj.__dict__.get(self.name, _SENTINEL)
        if value is not _SENTINEL:
            return value
        value = _freeze(self.func(obj))
        obj.__dict__[self.name] = value
        return value


def cached_method(func):
    """Per-instance cache for no-arg methods. Cleared by ``flush_cache()``."""
    name = f'__cached_method_{func.__name__}'

    @wraps(func)
    def wrapper(self):
        value = self.__dict__.get(name, _SENTINEL)
        if value is not _SENTINEL:
            return value
        value = _freeze(func(self))
        self.__dict__[name] = value
        return value
    return wrapper


def cached_args_method(func):
    """Per-instance, per-args-tuple cache. Args must be hashable."""
    name = f'__cached_args_method_{func.__name__}'

    @wraps(func)
    def wrapper(self, *args):
        cache = self.__dict__.get(name)
        if cache is None:
            cache = {}
            self.__dict__[name] = cache
        value = cache.get(args, _SENTINEL)
        if value is not _SENTINEL:
            return value
        value = _freeze(func(self, *args))
        cache[args] = value
        return value
    return wrapper


class class_cached_property:
    """Per-class cached property. Subclasses get their own entry.
    Host class must declare ``__class_cache__ = {}``.

    Slot-safe via ``getattr(obj, '__dict__', None)`` — this is the path
    ``cachedmethods 0.2.0`` broke (priority-rules bug #5)."""

    def __init__(self, func):
        self.__doc__ = getattr(func, '__doc__')
        self.func = func
        name = func.__name__
        if name.startswith('__') and not name.endswith('__'):
            name = f'_{func.__qualname__.split(".")[-2]}{name}'
        self.name = name

    def __get__(self, obj, cls):
        if obj is None:
            return self
        obj_dict = getattr(obj, '__dict__', None)
        if obj_dict is not None:
            value = obj_dict.get(self.name, _SENTINEL)
            if value is not _SENTINEL:
                return value
        class_cache = cls.__class_cache__.get(cls)
        if class_cache is None:
            class_cache = {}
            cls.__class_cache__[cls] = class_cache
        value = class_cache.get(self.name, _SENTINEL)
        if value is _SENTINEL:
            value = _freeze(self.func(obj))
            class_cache[self.name] = value
        if obj_dict is not None:
            obj_dict[self.name] = value
        return value


__all__ = [
    'FrozenDict',
    'cached_args_method',
    'cached_method',
    'cached_property',
    'class_cached_property',
]
