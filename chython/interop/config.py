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
"""
Configuration for the external tools.  `clean2d_engine` is `depict`'s, in `chython.depict._config`.

    >>> from chython.interop import config
    >>> config.conformer_engine = 'cdpkit'
    >>> config.class_paths = ['/opt/cdk.jar', '/opt/opsin.jar']
"""
from os import getenv
from typing import Literal


#: Conformer generation engine.  See `chython.interop.conformers.generate_conformers`.
conformer_engine: Literal['rdkit', 'cdpkit'] = 'rdkit'

#: JVM classpath for the Java tools, or `None` to take it from `CDK_PATH` and `OPSIN_PATH`.
#: `None` rather than a list built at import time, so the environment is read when the JVM starts and
#: not when `chython` is first imported.
class_paths: list[str] | None = None

_CONFORMER_ENGINES = ('rdkit', 'cdpkit')


class _Config(type(__import__('sys').modules[__name__])):
    """
    The module's own type, so that a bad assignment fails at the assignment.

    PEP 562 gives a module `__getattr__` but no `__setattr__`; replacing `__class__` on the module
    object gives both.  Used for one thing: refusing an engine name no converter implements, at the
    line that wrote it rather than inside a conformer call later.
    """
    __slots__ = ()

    def __getattribute__(self, name):
        if name == 'class_paths':
            explicit = super().__getattribute__('class_paths')
            if explicit is not None:
                return explicit
            return [getenv('CDK_PATH', 'cdk.jar'), getenv('OPSIN_PATH', 'opsin.jar')]
        return super().__getattribute__(name)

    def __setattr__(self, name, value):
        if name == 'conformer_engine' and value not in _CONFORMER_ENGINES:
            raise ValueError(f'conformer_engine must be one of {_CONFORMER_ENGINES}, got {value!r}')
        super().__setattr__(name, value)


__import__('sys').modules[__name__].__class__ = _Config


def _facade_alias(module_name: str, /, *names: str):
    """
    Make `names` on `module_name` live aliases of the same names here.

    An alias and not a copy: a copy reads correctly but breaks the write, since `chython.x = v` would
    rebind the facade's own name while every reader went on reading this module -- an assignment that
    appears to work and does nothing.
    """
    from sys import modules

    facade = modules[module_name]
    aliased = frozenset(names)
    cls = type(facade)
    # A fresh subclass per facade, so two callers cannot fight over one `__getattr__`.
    facade.__class__ = type(f'_Aliased{cls.__name__}', (cls,), {
        '__getattr__': lambda _, name: getattr(modules[__name__], name) if name in aliased
        else (_ for _ in ()).throw(AttributeError(f'module {module_name!r} has no attribute {name!r}')),
        '__setattr__': lambda self, name, value: setattr(modules[__name__], name, value)
        if name in aliased else cls.__setattr__(self, name, value),
    })
    # The names must NOT exist in the facade's own dict, or `__getattr__` is never consulted.
    for name in names:
        facade.__dict__.pop(name, None)


__all__ = ['conformer_engine', 'class_paths']
