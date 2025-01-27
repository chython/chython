# -*- coding: utf-8 -*-
#
#  Copyright 2024 Denis Lipatov <denis.lipatov163@gmail.com>
#  Copyright 2024 Vyacheslav Grigorev <slavick2000@yandex.ru>
#  Copyright 2024 Timur Gimadiev <timur.gimadiev@gmail.com>
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
This module introduces the `Vector` and `Polygon` classes, designed to perform mathematical 
calculations relevant to two-dimensional Cartesian coordinate systems and regular polygons. 

The `Vector` class facilitates operations with coordinates, including vector arithmetic 
(addition, subtraction, multiplication/division by scalars), normalization, rotation, and 
distance calculations. It also supports methods for determining the angle of a vector, its 
length, and whether it lies in a certain quadrant. Additionally, it includes functions for 
reflecting vectors about lines, finding the closest atom or point, and rotating vectors around 
other vectors or points.

The `Polygon` class focuses on properties and calculations related to regular polygons, such as 
finding the circumradius, calculating central angles, determining the apothem, and identifying 
the type of polygon based on its number of sides. It also includes static methods for 
calculating normals to lines defined by two vectors and for adding, averaging, and mirroring 
vectors.

Together, these classes provide a comprehensive toolkit for performing geometric computations 
essential in fields such as computer graphics, physics simulations, and computational chemistry, 
where precise manipulation and analysis of spatial relationships are required.
"""
import math
from typing import Union, TYPE_CHECKING, List

if TYPE_CHECKING:
    from .Properties import AtomProperties


class Vector:
    """
    The `Vector` class facilitates operations with coordinates, including vector arithmetic 
    (addition, subtraction, multiplication/division by scalars), normalization, rotation, and 
    distance calculations. It also supports methods for determining the angle of a vector, its 
    length, and whether it lies in a certain quadrant. Additionally, it includes functions for 
    reflecting vectors about lines, finding the closest atom or point, and rotating vectors around 
    other vectors or points.
    """
    def __init__(self, x: Union[int, float], y: Union[int, float]) -> None:
        """
        Constructor of the vector class

        Parameters:
        :param x Union[int, float]:
            The coordinate of the vector along the abscissa axis
        :param y Union[int, float]: 
            The coordinate of the vector along the ordinate axis

        Attributes:       
        x: the coordinate of the vector along the abscissa axis
        y: the coordinate of the vector along the ordinate axis
        """
        self.x: float = float(x)
        self.y: float = float(y)


    def __repr__(self) -> str:
        """
        The method needed for debugging the code

        Returns a string containing the coordinates of the point
        """
        return str(self.x) + ', ' + str(self.y)


    def copy(self) -> 'Vector':
        """
        Creates a copy of the current class object

        Returns a copy of the object
        """
        return Vector(self.x, self.y)


    def subtract(self, vector: 'Vector'):
        """
        A method for the operation of subtraction between vectors

        Parameters:
        :param vector 'Vector':
            Another object of the current class
        """
        self.x -= vector.x
        self.y -= vector.y


    def rotate(self, angle: float) -> None:
        """
        A method that rotates the vector by the appropriate angle from the signature 
        of the function and updates the coordinates of the current class object

        Parameters:
        :param angle float:
            The angle by which the vector should be rotated
        """
        new_x: float = self.x * math.cos(angle) - self.y * math.sin(angle)
        new_y: float = self.x * math.sin(angle) + self.y * math.cos(angle)

        self.x = new_x
        self.y = new_y


    def add(self, vector: 'Vector') -> None:
        """
        A class method that adds vectors and updates the coordinates of the current class object

        Parameters:
        :param vector 'Vector':
            Another object of the current class
        """
        self.x += vector.x
        self.y += vector.y


    def invert(self) -> None:
        """
        A class method that inverts the current coordinates of objects of the class
        """
        self.x = self.x * -1
        self.y = self.y * -1


    def divide(self, scalar: float) -> None:
        """
        A class method that divides the coordinates of the current class object
        vectors for an arbitrary number
        
        Parameters: 
        :param scalar float:
            Number divider
        """
        self.x = self.x / scalar
        self.y = self.y / scalar


    def normalise(self) -> None:
        """
        Normalization of coordinates (dividing them by the length of the vector itself)
        """
        if self.length() != 0:
            self.divide(self.length())


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
        return math.sqrt((self.x**2) + (self.y**2))


    def multiply_by_scalar(self, scalar: float) -> None:
        """
        Multiplies the coordinates of the current vector by an arbitrary real number
        
        Parameters: 
        :param scalar float
        """
        self.x = self.x * scalar
        self.y = self.y * scalar


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
        self.x -= vector.x
        self.y -= vector.y

        x = self.x * math.cos(angle) - self.y * math.sin(angle)
        y = self.x * math.sin(angle) + self.y * math.cos(angle)

        self.x = x + vector.x
        self.y = y + vector.y


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
        return (vector.x - self.x) ** 2 + (vector.y - self.y) ** 2


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
        a: float = self.y * vector.x
        b: float = self.x * vector.y

        if a > b:
            return 'clockwise'
        elif a == b:
            return 'neutral'
        else:
            return 'counterclockwise'


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

        a = (dx * dx - dy * dy) / (dx * dx + dy * dy)
        b = 2 * dx * dy / (dx * dx + dy * dy)

        new_x = a * (self.x - line_point_1.x) + b * (self.y - line_point_1.y) + line_point_1.x
        new_y = b * (self.x - line_point_1.x) - a * (self.y - line_point_1.y) + line_point_1.y

        self.x = new_x
        self.y = new_y


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
        d = (vector.x - vector_start.x) * (vector_end.y - vector_start.y) - (vector.y - vector_start.y) * (vector_end.x - vector_start.x)
        if d > 0:
            return 1
        elif d < 0:
            return -1
        else:
            return 0


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
        if determinant < 0:
            return 'clockwise'
        elif determinant == 0:
            return None
        else:
            return 'counterclockwise'


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
        difference = Vector.subtract_vectors(point_2, point_1)
        return difference.angle()


    @staticmethod
    def subtract_vectors(vector_1: 'Vector', vector_2: 'Vector')-> 'Vector':
        """
        Subtracts one vector from another.

        Parameters:
            vector_1 'Vector': 
                The vector from which to subtract.
            vector_2 'Vector': 
                The vector to subtract.

        Returns Vector: 
            A new Vector representing the result of the subtraction (vector_1 - vector_2).
        """
        x = vector_1.x - vector_2.x
        y = vector_1.y - vector_2.y
        return Vector(x, y)


    @staticmethod
    def add_vectors(vector_1: 'Vector', vector_2: 'Vector') -> 'Vector':
        """
        Adds two vectors together.

        Parameters:
            vector_1 'Vector': 
                The first vector to add.
            vector_2 'Vector': 
                The second vector to add.

        Returns Vector: 
            A new Vector representing the result of the addition (vector_1 + vector_2).
        """
        x = vector_1.x + vector_2.x
        y = vector_1.y + vector_2.y
        return Vector(x, y)


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
        x = (vector_1.x + vector_2.x) / 2
        y = (vector_1.y + vector_2.y) / 2
        return Vector(x, y)


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
        average_x = 0.0
        average_y = 0.0
        for vector in vectors:
            average_x += vector.x
            average_y += vector.y
        return Vector(average_x / len(vectors), average_y / len(vectors))


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
        delta = Vector.subtract_vectors(vector_2, vector_1)
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
        v1_x_diff: float = vector_1.x - origin.x
        v1_y_diff: float = vector_1.y - origin.y
        v2_x_diff: float = vector_2.x - origin.x
        v2_y_diff: float = vector_2.y - origin.y

        dot_product: float = v1_x_diff * v2_x_diff + v1_y_diff * v2_y_diff
        length_v1: float = math.sqrt(v1_x_diff ** 2 + v1_y_diff ** 2)
        length_v2: float = math.sqrt(v2_x_diff ** 2 + v2_y_diff ** 2)

        cos_angle = dot_product / (length_v1 * length_v2)
        return math.acos(cos_angle)













class Polygon:
    """
    The `Polygon` class focuses on properties and calculations related to regular polygons, such as 
    finding the circumradius, calculating central angles, determining the apothem, and identifying 
    the type of polygon based on its number of sides. It also includes static methods for 
    calculating normals to lines defined by two vectors and for adding, averaging, and mirroring 
    vectors.
    """
    def __init__(self, edge_number: float) -> None:
        """
        Initializes a Polygon with a specified number of edges.

        Parameters:
        :param edge_number float: 
                The number of edges (sides) of the polygon.
        """
        self.edge_number: float = edge_number

    @staticmethod
    def find_polygon_radius(edge_length: float, edge_number: float) -> float:
        """
        Calculates the radius of the circumcircle of a regular polygon.

        Parameters:
        :param edge_length float: 
            The length of one edge of the polygon.
        :param edge_number float: 
            The number of edges (sides) of the polygon.

        Returns float: 
            The radius of the circumcircle.
        """
        return edge_length / (2 * math.sin(math.pi / edge_number))

    @staticmethod
    def get_central_angle(edge_number: float) -> float:
        """
        Calculates the central angle of a regular polygon.

        Parameters:
        :param edge_number float: 
            The number of edges (sides) of the polygon.

        Returns float:
            The central angle in radians.
        """
        return math.radians(float(360) / edge_number)

    @staticmethod
    def get_apothem(radius: float, edge_number: float) -> float:
        """
        Calculates the apothem of a regular polygon.

        Parameters:
        :param radius float: 
            The radius of the circumcircle of the polygon.
        :param edge_number float: 
            The number of edges (sides) of the polygon.

        Returns float:
            The length of the apothem.
        """
        return radius * math.cos(math.pi / edge_number)

    @staticmethod
    def get_apothem_from_side_length(length: float, edge_number: float) -> float:
        """
        Calculates the apothem of a regular polygon given the side length.

        Parameters:
        :param length float: 
            The length of one edge of the polygon.
        :param edge_number float: 
            The number of edges (sides) of the polygon.

        Returns float: 
            The length of the apothem.
        """
        radius: float = Polygon.find_polygon_radius(length, edge_number)
        return Polygon.get_apothem(radius, edge_number)