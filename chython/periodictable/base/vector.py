# -*- coding: utf-8 -*-
#
#  Copyright 2024 Denis Lipatov <denis.lipatov163@gmail.com>
#  Copyright 2024 Vyacheslav Grigorev <slavick2000@yandex.ru>
#  Copyright 2024 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2024, 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from math import cos, sin, hypot, atan2


class Vector:
    __slots__ = ('x', 'y')

    def __init__(self, x: float = 0., y: float = 0.):
        self.x = x
        self.y = y

    def __repr__(self):
        return f'Vector({self.x}, {self.y})'

    def __neg__(self):
        """
        A class method that inverts the current coordinates of objects of the class
        """
        return Vector(-self.x, -self.y)

    def __sub__(self, vector: 'Vector'):
        """
        A method for the operation of subtraction between vectors
        """
        return Vector(self.x - vector.x, self.y - vector.y)

    def __add__(self, vector: 'Vector'):
        """
        A method for the operation of addition between vectors
        """
        return Vector(self.x + vector.x, self.y + vector.y)

    def __truediv__(self, scalar: float):
        """
        A class method that divides the coordinates of the vector by a given scalar
        """
        return Vector(self.x / scalar, self.y / scalar)

    def __mul__(self, scalar: float):
        """
        Multiplies the coordinates of the current vector by an arbitrary real number
        """
        return Vector(self.x * scalar, self.y * scalar)

    def __float__(self):
        """
        Calculates the length of the current vector

        Returns float
        """
        return hypot(self.x, self.y)

    def __iter__(self):
        yield self.x
        yield self.y

    def __len__(self):
        return 2

    def __matmul__(self, vector: 'Vector'):
        return self.x * vector.y - self.y * vector.x

    def __or__(self, vector: 'Vector'):
        """
        Calculate distance between two vectors
        """
        return hypot(vector.x - self.x, vector.y - self.y)

    def rotate(self, angle: float, vector: 'Vector' = None):
        """
        A method that rotates the vector by the angle in radians
        """
        c = cos(angle)
        s = sin(angle)
        if vector is None:
            return Vector(self.x * c - self.y * s, self.x * s + self.y * c)
        xy = self - vector
        return vector + Vector(xy.x * c - xy.y * s, xy.x * s + xy.y * c)

    def normalise(self):
        """
        Normalization of coordinates (dividing them by the length of the vector itself)
        """
        if ln := float(self):
            return self / ln
        return self

    def angle(self, vector: 'Vector' = None) -> float:
        """
        A method calculates the angle of inclination of the current vector
        or the vector between given vector and the current vector.
        """
        if vector is None:
            return atan2(self.y, self.x)
        else:
            return atan2(vector.y - self.y, vector.x - self.x)


__all__ = ['Vector']
