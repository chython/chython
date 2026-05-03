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
from functools import cache


@cache
def _start_jvm():
    """Start JVM once with all Java dependencies on classpath. Thread-safe via @cache."""
    from jpype import isJVMStarted, startJVM

    if not isJVMStarted():
        from chython import class_paths

        startJVM('--enable-native-access=ALL-UNNAMED', classpath=class_paths)


@cache
def get_cdk():
    """Get CDK Java package. Starts JVM if needed."""
    try:
        from jpype import JPackage

        _start_jvm()
        return JPackage('org').openscience.cdk
    except (ImportError, AttributeError):
        raise ImportError('Java/JPype/CDK.jar is not installed or broken. make sure CDK_PATH env variable is set')


@cache
def get_opsin():
    """Get OPSIN NameToStructure instance. Starts JVM if needed."""
    try:
        from jpype import JPackage

        _start_jvm()
        return JPackage('uk').ac.cam.ch.wwmm.opsin.NameToStructure.getInstance()
    except (ImportError, AttributeError):
        raise ImportError('Java/JPype/OPSIN.jar is not installed or broken. make sure OPSIN_PATH env variable is set')


__all__ = ['get_cdk', 'get_opsin']
