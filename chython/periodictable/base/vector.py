# -*- coding: utf-8 -*-
#
#  Copyright 2024, 2025 Denis Lipatov <denis.lipatov163@gmail.com>
#  Copyright 2024, 2025 Vyacheslav Grigorev <slavick2000@yandex.ru>
#  Copyright 2024, 2025 Timur Gimadiev <timur.gimadiev@gmail.com>
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

"""
This module introduces the `Vector` class, designed to perform mathematical 
calculations relevant to two-dimensional Cartesian coordinate systems. 

The `Vector` class facilitates operations with coordinates, including vector arithmetic 
(addition, subtraction, multiplication/division by scalars), normalization, rotation, and 
distance calculations. It also supports methods for determining the angle of a vector, its 
length, and whether it lies in a certain quadrant. Additionally, it includes functions for 
reflecting vectors about lines, finding the closest atom or point, and rotating vectors around 
other vectors or points.

This class provides a comprehensive toolkit for performing geometric computations 
essential in fields such as computer graphics, physics simulations, and computational chemistry, 
where precise manipulation and analysis of spatial relationships are required.
"""
import numpy as np
import math
from typing import List, TYPE_CHECKING

if TYPE_CHECKING:
    from ...algorithms.calculate2d.Properties import *

class Vector(np.ndarray):
    """
    The `Vector` class facilitates operations with coordinates, including vector arithmetic 
    (addition, subtraction, multiplication/division by scalars), normalization, rotation, and 
    distance calculations. It also supports methods for determining the angle of a vector, its 
    length, and whether it lies in a certain quadrant. Additionally, it includes functions for 
    reflecting vectors about lines, finding the closest atom or point, and rotating vectors around 
    other vectors or points.
    """
    __slots__ = ()


    def __new__(cls, *args) -> np.ndarray:
        """
        Constructor of the Vector class.

        Parameters:
        :param args: A variable number of arguments representing the coordinates of the vector.
                    Each argument in args is expected to be a number (Union[int, float]).
                    The first two arguments correspond to the x and y coordinates, respectively.

        Returns:
        Vector: An instance of the Vector class initialized with the provided coordinates.
        """
        obj = np.asarray(args, dtype=np.float64).view(cls)
        return obj


    @property
    def x(self):
        """
        Gets the x-coordinate of the vector along the abscissa axis.
        :return: The x-coordinate value (type: Union[int, float]).
        """
        return self[0]
    
    
    @property
    def y(self):
        """
        Gets the y-coordinate of the vector along the ordinate axis.
        :return: The y-coordinate value (type: Union[int, float]).
        """
        return self[1]
    
  
    def __repr__(self):
        return f'Vector({self.x}, {self.y})'


    def __float__(self):
        """
        Calculates the length of the current vector

        Returns float
        """
        return math.hypot(self.x, self.y)


    def __iter__(self):
        yield self.x
        yield self.y


    def __len__(self):
        return 2


    def __or__(self, vector: 'Vector'):
        """
        Calculate distance between two vectors
        """
        return math.hypot(vector.x - self.x, vector.y - self.y)
    
    
    def rotate(self, angle: float) -> None:
        """
        A method that rotates the vector by the appropriate angle from the signature 
        of the function and updates the coordinates of the current class object

        Parameters:
        :param angle float:
            The angle by which the vector should be rotated
        """
        rotation_matrix = np.array([
            [np.cos(angle), -np.sin(angle)],
            [np.sin(angle), np.cos(angle)]
        ], dtype=np.float64)
        self[:] = rotation_matrix @ self


    def invert(self) -> None:
        """
        A class method that inverts the current coordinates of objects of the class
        """
        self[:] = -self


    def normalise(self) -> None:
        """
        Normalization of coordinates (dividing them by the length of the vector itself)
        """
        if length := np.linalg.norm(self):
            self[:] /= length


    def angle(self) -> float:
        """
        A method that calculates the angle of inclination of the current vector

        Returns float the angle of inclination of the vector
        """
        return math.atan2(self.y, self.x)


    def length(self) -> float:
        """
        Calculates the length of the current vector
        
        Returns float 
        """
        return np.linalg.norm(self)


    def rotate_around_vector(self, angle: float, vector: 'Vector') -> None:
        """
        Rotates a point (or vector) around a given vector by a specified angle.

        Parameters:
        :param angle float: 
            The angle by which to rotate the point, typically measured in radians.
        :param vector 'Vector':
            The vector around which the rotation occurs. This vector serves as the reference 
            point.
        """
        translated_vector = self - vector
        rotation_matrix = np.array([
            [np.cos(angle), -np.sin(angle)],
            [np.sin(angle), np.cos(angle)]
        ])
        rotated_vector = rotation_matrix @ translated_vector
        self[:] = np.add(rotated_vector, vector)


    def get_closest_atom(self, atom_1: 'AtomProperties', atom_2: 'AtomProperties') -> 'AtomProperties':
        """
        This method determines which of the two atoms (represented by the objects atom_1 and atom_2)
        is closer to the current object (represented by self).
        
        Parameters:
        :param atom_1: 'AtomProperties':
            The first atom to compare.
        :param atom_2: 'AtomProperties': 
            The second atom to compare.
        
        Returns 'AtomProperties':
            The closest atom.
        """
        distance_1 = self.get_squared_distance(atom_1.position)
        distance_2 = self.get_squared_distance(atom_2.position)
        return atom_1 if distance_1 < distance_2 else atom_2


    def get_closest_point_index(self, point_1: 'Vector', point_2: 'Vector') -> int:
        """
        The method is designed to determine which of the two specified coordinates (point_1: 'Vector', point_2: 'Vector')
        closer to the current point.

        Parameters
        :param point_1 'Vector': 
            The first point to be compared with. It can be a tuple, a list, or an object 
            representing coordinates.
        :param point_2 'Vector': 
            The second point to compare with. Similarly, it can be a tuple, a list, or an 
            object.        

        Returns int:
            The index of the nearest point: 0 for point_1 and 1 for point_2.
        """
        distance_1 = self.get_squared_distance(point_1)
        distance_2 = self.get_squared_distance(point_2)
        return 0 if distance_1 < distance_2 else 1


    def get_squared_length(self) -> float:
        """
        Calculates the length squared

        Returns float:
            Vector length squared
        """
        return self.x ** 2 + self.y ** 2


    def get_squared_distance(self, vector: 'Vector') -> float:
        """
        The method is designed to calculate the square of the distance between the current vector
        (represented by self) and the specified vector (or point) represented by the vector object. 
        
        Parameters
        :param vector: 'Vector': 
            An object representing a vector or point from which to calculate the distance.
        
        Returns float:
            The square of the distance
        """
        return np.sum((self - vector) ** 2)


    def get_distance(self, vector: 'Vector') -> float:
        """
        The method is designed to calculate the distance between the current vector (represented by self) and 
        the specified vector (or point) represented by the vector object. 
        
        Parameters
        :param vector: 'Vector': 
            An object representing a vector or point from which to calculate the distance.
        
        Returns float:
            The distance between the coordinates of the current vector and the passed parameter
        """
        return math.sqrt(self.get_squared_distance(vector))


    def get_rotation_away_from_vector(self, vector: 'Vector', center: 'Vector', angle: float) -> float:
        """
        The method is designed to determine how much the angle of rotation (in a positive or negative direction)
        from a given vector measures the distance to this vector. 
        
        Parameters
        :param vector 'Vector': 
            The vector to "move away from". It can be a point or a direction, relative to which 
            the rotation is taking place.
        :param center 'Vector': 
            The center of rotation around which the object (represented by self) rotates.
        :param angle float: 
            The angle at which the rotation occurs. This value can be positive or negative.
        
        Returns returns the rotation angle that minimizes the distance to the vector,
            either in a positive or negative direction.
        """
        tmp = self.copy()
        tmp.rotate_around_vector(angle, center)
        squared_distance_1 = tmp.get_squared_distance(vector)
        tmp.rotate_around_vector(-2.0 * angle, center)
        squared_distance_2 = tmp.get_squared_distance(vector)
        return angle if squared_distance_2 < squared_distance_1 else -angle


    def rotate_away_from_vector(self, vector: 'Vector', center: 'Vector', angle: float) -> None:
        """
        The method is designed to rotate the current object (represented by self) around a given 
        one center in such a way as to minimize the distance to the specified vector. 
        If rotation in one direction leads to a decrease in the distance, the function corrects 
        the rotation,to ensure maximum distance from the vector.

        Parameters
        :param vector 'Vector': 
            The vector to "move away from". It can be a point or a direction, relative to which 
            the rotation is taking place.
        :param center 'Vector': 
            The center of rotation around which the object rotates. 
        :param angle float:
            The angle at which the rotation occurs. This value can be positive or negative.
        """
        self.rotate_around_vector(angle, center)
        squared_distance_1 = self.get_squared_distance(vector)
        self.rotate_around_vector(-2.0 * angle, center)
        squared_distance_2 = self.get_squared_distance(vector)
        if squared_distance_2 < squared_distance_1:
            self.rotate_around_vector(2.0 * angle, center)


    def get_clockwise_orientation(self, vector: 'Vector') -> str:
        """
        The method is designed to determine the orientation (positive or negative) between
        the current object (represented by self) and the specified vector (represented by 
        the vector object). 
        
        Parameters
        :param vector 'Vector':
            The vector relative to which the orientation is determined.

        Returns str:
            A string indicating whether the orientation is "clockwise", "counterclockwise" 
            or "neutral".
        """
        orientation_value = self.y * vector.x - self.x * vector.y
        return 'clockwise' if orientation_value > 0 else 'counterclockwise' \
            if orientation_value < 0 else 'neutral'


    def mirror_about_line(self, line_point_1: 'Vector', line_point_2: 'Vector') -> None:
        """
        The method is designed to reflect the current object (represented by self) relative to a 
        given line, defined by two points (line_point_1 and line_point_2). After performing this 
        function, the coordinates of the object will be changed so that it is on the opposite 
        side of the line, keeping the same distance to the line.

        Parameters
        :param line_point_1: 'Vector':
            The first point defining the line. 
        :param line_point_2: 'Vector': 
            The second point defining the line.
        """
        dx = line_point_2.x - line_point_1.x
        dy = line_point_2.y - line_point_1.y

        a = (dx**2 - dy**2) / (dx**2 + dy**2)
        b = 2 * dx * dy / (dx**2 + dy**2)
        self[:] = np.array([
            a * (self.x - line_point_1.x) + b * (self.y - line_point_1.y) + line_point_1.x,
            b * (self.x - line_point_1.x) - a * (self.y - line_point_1.y) + line_point_1.y
        ])


    @staticmethod
    def get_position_relative_to_line(vector_start: 'Vector', vector_end: 'Vector', vector: 'Vector') -> int:
        """
        Determines the position of a vector relative to a line defined by two points.
        
        Parameters:
        :param vector_start 'Vector':
            The start point of the line.
        :param vector_end 'Vector': 
            The end point of the line.
        :param vector 'Vector': 
            The vector whose position relative to the line is to be determined.
        
        Returns int:
            1 if the vector is to the left of the line, -1 if the vector is to the right of the 
            line, 0 if the vector lies on the line.
        """
        determinant = (vector.x - vector_start.x) * (vector_end.y - vector_start.y) - \
                    (vector.y - vector_start.y) * (vector_end.x - vector_start.x)
        return np.sign(determinant)


    @staticmethod
    def get_directionality_triangle(vector_a: 'Vector', vector_b: 'Vector', vector_c: 'Vector') -> str:
        """
        Determines the directionality of the triangle formed by three vectors (or points).

        Parameters:
        :param vector_a 'Vector': 
            The first vertex of the triangle.
        :param vector_b 'Vector': 
            The second vertex of the triangle.
        :param vector_c 'Vector': 
            The third vertex of the triangle.

        Returns str: 
                - 'clockwise' if the triangle is oriented in a clockwise direction.
                - 'counterclockwise' if the triangle is oriented in a counterclockwise direction.
                - None if the three points are collinear (lie on the same line).
        """
        determinant = (vector_b.x - vector_a.x) * (vector_c.y - vector_a.y) - \
                      (vector_c.x - vector_a.x) * (vector_b.y - vector_a.y)
        if determinant:
            return 'clockwise'  if determinant < 0  else 'counterclockwise'


    @staticmethod
    def mirror_vector_about_line(line_point_1: 'Vector', line_point_2: 'Vector', point: 'Vector')-> 'Vector':
        """
        Mirrors a point (or vector) across a line defined by two points.

        Parameters:
        :param line_point_1 'Vector': 
            The first point defining the line.
        :param line_point_2 'Vector': 
            The second point defining the line.
        :param point 'Vector': 
            The point to be mirrored across the line.

        Returns Vector:
            A new Vector representing the mirrored point across the line.
        """
        dx = line_point_2.x - line_point_1.x
        dy = line_point_2.y - line_point_1.y
        a = (dx * dx - dy * dy) / (dx * dx + dy * dy)
        b = 2 * dx * dy / (dx * dx + dy * dy)
        x_new = a * (point.x - line_point_1.x) + b * (point.y - line_point_1.y) + line_point_1.x
        y_new = b * (point.x - line_point_1.x) - a * (point.y - line_point_1.y) + line_point_1.y
        return Vector(x_new, y_new)


    @staticmethod
    def get_line_angle(point_1: 'Vector', point_2: 'Vector') -> float:
        """
        Calculates the angle of a line defined by two points with respect to the positive x-axis.

        Parameters:
            point_1 'Vector': 
                The first point defining the line.
            point_2 'Vector':
                The second point defining the line.

        Returns float:
            The angle of the line in radians, in the range [-π, π].
        """
        difference: 'Vector' = point_2 - point_1
        return difference.angle()


    @staticmethod
    def get_midpoint(vector_1: 'Vector', vector_2: 'Vector') -> 'Vector':
        """
        Calculates the midpoint between two vectors.

        Parameters:
            vector_1 'Vector': 
                The first vector.
            vector_2 'Vector': 
                The second vector.

        Returns Vector: 
            A new Vector representing the midpoint between vector_1 and vector_2.
        """
        midpoint: Vector = (vector_1 + vector_2) / 2
        return midpoint


    @staticmethod
    def get_average(vectors: List['Vector']) -> 'Vector':
        """
        Calculates the average of a list of vectors.

        Parameters:
        :param vectors List[Vector]:
            A list of vectors for which the average is to be calculated.

        Returns:
            Vector: A new Vector representing the average of the input vectors.
        """
        vectors_array = np.array(vectors)  
        average = np.mean(vectors_array, axis=0)
        return Vector(*average)  


    @staticmethod
    def get_normals(vector_1: 'Vector', vector_2: 'Vector') -> List['Vector']:
        """
        Calculates the normal vectors to the line defined by two vectors.

        Parameters:
        :param vector_1 'Vector': 
                The first vector defining the line.
        :param vector_2 'Vector': 
                The second vector defining the line.

        Returns List[Vector]: 
            A list containing two normal vectors to the line defined by vector_1 and vector_2.
        """
        delta: Vector = vector_2 - vector_1
        return [Vector(-delta.y, delta.x), Vector(delta.y, -delta.x)]


    @staticmethod
    def get_angle_between_vectors(vector_1: 'Vector', vector_2: 'Vector', origin: 'Vector') -> float:
        """
        Calculates the angle between two vectors relative to a given origin point.

        Parameters:
        :param vector_1 'Vector': 
            The first vector.
        :param vector_2 'Vector': 
            The second vector.
        :param origin 'Vector': 
            The origin point relative to which the angle is calculated.

        Returns:
            float: The angle between vector_1 and vector_2 in radians, in the range [0, π].
        """
        v1_diff: np.array = np.subtract(vector_1, origin)
        v2_diff: np.array = np.subtract(vector_2, origin)
        dot_product = np.dot(v1_diff, v2_diff)

        length_v1: float = np.linalg.norm(v1_diff)
        length_v2: float = np.linalg.norm(v2_diff)

        cos_angle: float = dot_product / (length_v1 * length_v2)
        cos_angle: float = np.clip(cos_angle, -1.0, 1.0)
        return np.arccos(cos_angle)
        
__all__ = ['Vector']