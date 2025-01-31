# -*- coding: utf-8 -*-
#
#  Copyright 2024, 2025 Denis Lipatov <denis.lipatov163@gmail.com>
#  Copyright 2024, 2025 Vyacheslav Grigorev <slavick2000@yandex.ru>
#  Copyright 2024, 2025 Timur Gimadiev <timur.gimadiev@gmail.com>
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
This module introduces the `Polygon` class, designed to perform mathematical 
calculations relevant to two-dimensional Cartesian coordinate systems and regular polygons. 

The `Polygon` class focuses on properties and calculations related to regular polygons, such as 
finding the circumradius, calculating central angles, determining the apothem, and identifying 
the type of polygon based on its number of sides. It also includes static methods for 
calculating normals to lines defined by two points and for adding, averaging, and mirroring points.

This class provides a comprehensive toolkit for performing geometric computations 
essential in fields such as computer graphics, physics simulations, and computational chemistry, 
where precise manipulation and analysis of spatial relationships is required.
"""
import math
    
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
    
__all__ = ['Polygon']
    