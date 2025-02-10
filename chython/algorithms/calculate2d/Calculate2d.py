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
Class for calculating the 2D layout of a molecular graph, returning the coordinates of atom 
vertices in a molecular container.

This class provides methods for calculating and optimizing the 2D structure of a molecule, 
including determining atom coordinates, handling collisions, and defining properties of rings 
and bonds between atoms. Key functions include:

- Calculating the initial positions of atoms and their subsequent adjustment to minimize 
overlaps.
- Defining and classifying rings within the molecule, including handling bridged, spiro-fused, 
and condensed rings.
- Handling cis-trans isomerism and atom configurations.
- Working with various types of bonds (single, double, triple) and their impact on atom 
orientation.
- Calculating atom positions in ring structures and aromatic compounds.
- Handling collisions between atoms to improve molecule visualization.

The class uses auxiliary functions for vector operations, determining atom neighbors, 
calculating angles, and handling specific cases such as atoms with one, two, three, or four 
neighbors. It also includes methods for working with ring structures, including determining the 
ring center and positioning atoms within the ring.

Utilizes AtomProperties, BondProperties, RingProperties, RingOverlap classes to represent atoms, 
bonds, and rings, respectively. Additionally, the KKLayout class, which implements the 
Kamada-Kawai algorithm, is used to minimize the system's energy by representing atoms as masses 
connected by springs with a certain stiffness (used for calculating coordinates in bridged 
cyclic molecules).
"""
from typing import List, Dict, Optional, Set, Tuple, Union, Generator, TYPE_CHECKING
from ...periodictable.base.vector import Vector
from .polygon import Polygon
from .KKLayout import KKLayout
from .Properties import * # RingProperties, AtomProperties, RingOverlap
import math # radiand, sin, cos, pi
import numpy as np

if TYPE_CHECKING:
    from ...containers import MoleculeContainer
    from ...containers import Bond


class Calculate2d:
    """
    Class for calculating the 2D layout of a molecular graph, returning the coordinates of atom 
    vertices in a molecular container.
    """

    def __init__(self) -> None:
        """
        The initial attributes initialization of the class includes:
        bond_length: int: The bond length between atoms in the molecular graph. Used to determine 
            the distance between atoms when calculating their coordinates.
        overlap_sensitivity: float: Sensitivity to atom overlap. Used to determine how close atoms 
            can be to each other before they are considered to overlap.
        overlap_resolution_iterations: int: The number of iterations for resolving atom overlaps.
            Indicates how many times the algorithm will attempt to improve atom positioning to 
            minimize overlaps.
        ring_overlaps: List['RingOverlap']: A list of RingOverlap objects representing overlaps 
            between rings in the molecule. Used for identifying and handling ring overlaps.
        total_overlap_score: float: The total overlap score in the molecule. Used for evaluating the 
            quality of atom positioning and minimizing overlaps.
        finetune: bool: A flag indicating the need for detailed adjustment (finetuning) of atom 
            positions after the main calculation. Currently always set to True, but can be changed 
            to disable detailed adjustment.
        ring_overlap_id_tracker: int: A counter for assigning unique identifiers to RingOverlap 
            objects. Used for the unique identification of ring overlaps.
        ring_id_tracker: int: A counter for assigning unique identifiers to rings. Simplifies ring 
            management in the structure.
        rings: List['RingProperties']: A list of RingProperties objects representing all rings in 
            the molecule. Used for working with ring structures and their attributes.
        id_to_ring: Dict[int, 'RingProperties']: A dictionary mapping ring identifiers to their 
            RingProperties objects. Facilitates access to ring properties by their identifiers.
        """
        
        self.bond_length: int = 15
        self.overlap_sensitivity: float = 0.10
        self.overlap_resolution_iterations: int = 5
        self.ring_overlaps: List['RingOverlap'] = []
        self.total_overlap_score: float = 0.0

        self.ring_overlap_id_tracker: int = 0
        self.ring_id_tracker: int = 0
        self.rings: List['RingProperties'] = []
        self.id_to_ring: Dict[int, 'RingProperties'] = {}

    def _calculate2d_coord(self, order: List[int],
                           mc: 'MoleculeContainer') -> List[List[float]]:
        self.create_property_attributes(mc)
        self.define_rings()

        self.initial_approximation()
        self.collision_handling()
        return self.get_coord(order)


    
    
    def create_property_attributes(self, mc: 'MoleculeContainer') -> None:
        """
        Initializes the property attributes for the molecular container, including atoms, bonds, and their 
        relationships.

        This method sets up the initial properties for the molecular container by creating dictionaries for 
        atoms and bonds based on the adjacency information provided by the molecular container. It also 
        refreshes the neighbours list for each atom to ensure accurate representation of the molecular 
        graph.

        Parameters
        :param mc: MoleculeContainer:
            An instance of MoleculeContainer that holds the molecular graph, including atoms, bonds, and 
            rings information necessary for the 2D layout calculation.

        Notes:
        - Initializes the `atoms` dictionary with atom indices as keys and AtomProperties instances as 
        values, where each AtomProperties instance is created with the atom index and its corresponding 
        symbol.
        - Refreshes the neighbours list for each atom based on the adjacency information from the molecular 
        container.
        - Initializes the `bonds` dictionary with tuples of atom indices as keys and BondProperties 
        instances as values, representing the bonds between atoms.
        - Sets up the `graph` dictionary to map each atom to its list of neighbouring atoms, facilitating 
        the representation of the molecular structure.

        The method is crucial for preparing the molecular container for further calculations by establishing 
        the basic properties and relationships between atoms and bonds, which are essential for the 2D 
        layout calculation.
        """
        self.mc: 'MoleculeContainer' = mc

        self.atoms: Dict[int, AtomProperties] = {}
        for atom in self.mc.int_adjacency:
            self.atoms[atom] = AtomProperties(atom)
        self.refresh_neighbours(self.mc.int_adjacency)

        self.graph: Dict['AtomProperties', List['AtomProperties']] = {}
        for atom in self.atoms.values():
            self.graph[atom] = atom.neighbours


    def refresh_neighbours(self, graph: Dict[int, Dict[int, int]]) -> None:
        """
        Refreshes the neighbours list for each atom based on the provided graph adjacency information.

        This method updates the neighbours list for each atom in the molecular graph to ensure an accurate 
        representation of the molecular structure. It iterates through the graph, which is a dictionary 
        mapping atom indices to their adjacent atom indices, and assigns the corresponding AtomProperties 
        instances to the neighbours list of each atom.

        Parameters
        :param graph: Dict[int, Dict[int, int]]:
            A dictionary where keys are atom indices and values are dictionaries mapping to adjacent atom 
            indices. This structure represents the adjacency information of the molecular graph, indicating
            which atoms are directly connected.
        """
        for atom_index, neighbor_indexes in graph.items():
            neighbours: List['AtomProperties'] = []
            for neighbour_index in neighbor_indexes:
                neighbours.append(self.atoms[neighbour_index])
            self.atoms[atom_index].neighbours = neighbours


    def bond_lookup(self, atom: 'AtomProperties', next_atom: 'AtomProperties') -> Optional['Bond']:
        """
        Retrieves the bond object between two specified atoms.

        This method checks if there is a bond between the given atoms and returns the corresponding 
        bond object if it exists.

        Parameters:
        atom (AtomProperties): The first atom involved in the bond.
        next_atom (AtomProperties): The second atom involved in the bond.

        Returns:
        Optional[Bond]: The bond object between the two atoms if it exists; otherwise, returns None.
        """
        if self.mc.has_bond(atom.id, next_atom.id):
            return self.mc.bond(atom.id, next_atom.id)

    def get_configuration(self, atom1: 'AtomProperties', atom2: 'AtomProperties') -> Optional[str]:
        """
        Checks if there is a configuration between the specified atoms.

        If no configuration exists, it returns None. Otherwise, it returns 
        'cis' if the configuration is cis or 'trans' if the configuration is trans.

        The configurations are stored in a dictionary, where the keys are tuples 
        of atom IDs that may have configurations. The condition is that the bond 
        between these atoms must be a double bond. The values in the dictionary 
        are boolean: True indicates a cis configuration, while False indicates 
        a trans configuration. If no configurations are defined, the dictionary 
        will be empty.

        Parameters:
        atom1 (AtomProperties): The first atom to check for configuration.
        atom2 (AtomProperties): The second atom to check for configuration.

        Returns:
        Optional[str]: 'cis' if the configuration is cis, 'trans' if it is trans,
                    or None if no configuration exists.
        """
        if (atom1.id, atom2.id) in list(self.mc._cis_trans_stereo.keys()):
            configuration = self.mc._cis_trans_stereo[(atom1.id, atom2.id)]
            return 'cis' if configuration else 'trans'

    def define_rings(self) -> None:
        """
        Defines the rings within the molecule, identifies ring overlaps, and handles bridged ring systems.

        This method performs several key steps in the process of analyzing the molecular structure:
        1. It initializes the rings present in the molecule by converting the simple cycle list (SSSR) from 
        the molecular container into RingProperties objects.
        2. Identifies overlaps between rings and creates RingOverlap objects for them.
        3. Finds and processes bridged ring systems, creating a unified representation for interconnected 
        rings that share atoms.

        Notes
        -----
        - Initially, it retrieves the simple cycle list (SSSR) from the molecular container and converts 
        each cycle into a RingProperties object, adding them to the class's ring list.
        - It then iterates through all pairs of rings to identify overlaps, creating RingOverlap objects for 
        those that share atoms and adding them to the class's ring overlaps list.
        - For each ring, it updates the list of neighbouring rings based on identified overlaps, enhancing 
        the representation of the molecular structure's connectivity.
        - The method also handles bridged ring systems by identifying rings that are part of a larger, 
        interconnected system and merges them into a single RingProperties object, ensuring a coherent 
        representation of complex cyclic structures.
        - This process involves finding all rings involved in a bridged system, removing the original rings 
        from the list, and adding a new RingProperties object that represents the bridged system.
        - Finally, it iterates through all rings to find any bridged systems not yet processed and repeats 
        the merging process, ensuring that all interconnected rings are represented as unified entities.
        """
        rings = self.mc.sssr
        if not rings:
            return
        
        for neighbor_indexes in rings:
            members_ring: List['AtomProperties'] = [self.atoms[atom_index] for atom_index in neighbor_indexes]
            ring = RingProperties(members_ring)
            self.add_ring(ring)

        for i, ring_1 in enumerate(self.rings[:-1]):
            for ring_2 in self.rings[i + 1:]:
                ring_overlap = RingOverlap(ring_1, ring_2)
                if len(ring_overlap.atoms) > 0:
                    self.add_ring_overlap(ring_overlap)

        for ring in self.rings:
            neighbouring_rings = self.find_neighbouring_rings(self.ring_overlaps, ring.id)
            ring.neighbouring_rings = neighbouring_rings

        while True:
            ring_id: int = -1
            for ring in self.rings:
                if self.is_part_of_bridged_ring(ring.id) and not ring.bridged:
                    ring_id: int = ring.id
            if ring_id == -1:
                break
            ring: 'RingProperties' = self.id_to_ring[ring_id]

            involved_ring_ids: Union[list[int], Set[int]] = []
            self.get_bridged_ring_subrings(ring.id, involved_ring_ids)
            involved_ring_ids = set(involved_ring_ids)

            self.has_bridged_ring = True
            self.create_bridged_ring(involved_ring_ids)

            for involved_ring_id in involved_ring_ids:
                involved_ring = self.id_to_ring[involved_ring_id]
                self.remove_ring(involved_ring)

        bridged_systems = self.find_bridged_systems(
            self.rings, self.ring_overlaps)
        if bridged_systems and not self.has_bridged_ring:
            self.has_bridged_ring = True
            for bridged_system in bridged_systems:
                involved_ring_ids = set(bridged_system)
                self.create_bridged_ring(involved_ring_ids)
                for involved_ring_id in involved_ring_ids:
                    involved_ring = self.id_to_ring[involved_ring_id]
                    self.remove_ring(involved_ring)


    
    def add_ring_overlap(self, ring_overlap: 'RingOverlap') -> None:
        """
        Adds a new ring overlap to the list of ring overlaps and assigns it a unique identifier.

        This method assigns a unique identifier to the given ring overlap and appends it to the class's list 
        of ring overlaps. It ensures that each ring overlap is uniquely identifiable and can be tracked 
        throughout the calculation process.

        Parameters
        :param ring_overlap: RingOverlap
            The ring overlap to be added to the list of overlaps. This object represents the intersection 
            between two rings in the molecule, which may need special handling during the layout calculation 
            to avoid visual clutter or incorrect representation.
        """
        ring_overlap.id = self.ring_overlap_id_tracker
        self.ring_overlaps.append(ring_overlap)
        self.ring_overlap_id_tracker += 1


    def is_part_of_bridged_ring(self, ring_id: int) -> bool:
        """
        Determines if a given ring is part of a bridged ring system.

        This method checks if a ring, identified by its ID, is involved in a bridged ring system by 
        examining the list of ring overlaps. It returns True if the ring is part of a bridged system, 
        indicating that it is connected to another ring through a bridge, and False otherwise.

        Parameters
        :param ring_id: int:
            The identifier of the ring to check for involvement in a bridged ring system.

        Returns bool:
            True if the ring is part of a bridged ring system, indicating that it is interconnected with 
            another ring through a bridge, and False otherwise.
        """
        return any(ring_overlap.involves_ring(ring_id) \
            and ring_overlap.is_bridge() for ring_overlap in self.ring_overlaps)


    def get_bridged_ring_subrings(self, ring_id: int, involved_ring_ids: List[int]) -> None:
        """
        Recursively identifies and collects the IDs of all rings involved in a bridged ring system starting 
        from a given ring ID.

        This method is used to find all rings that are interconnected as part of a bridged ring system, 
        starting from a specified ring ID. It recursively explores neighboring rings to identify all rings 
        that are connected through bridges, adding their IDs to a list of involved ring IDs.

        Parameters
        :param ring_id: int
            The identifier of the starting ring from which to begin the search for interconnected rings in a 
            bridged system.
        :param involved_ring_ids: List[int]
            A list to which the IDs of rings involved in the bridged system are appended. This list is 
            populated with the IDs of all rings found to be part of the bridged ring system.
        """
        involved_ring_ids.append(ring_id)
        ring = self.id_to_ring[ring_id]
        for neighbour_id in ring.neighbouring_rings:
            is_connected_by_bridge: bool = \
                neighbour_id not in involved_ring_ids and neighbour_id != ring_id and \
                self.rings_connected_by_bridge(self.ring_overlaps, ring_id, neighbour_id)
            if is_connected_by_bridge:
                self.get_bridged_ring_subrings(neighbour_id, involved_ring_ids)

    @staticmethod
    def rings_connected_by_bridge(ring_overlaps: List['RingOverlap'],
                                  ring_id_1: int, ring_id_2: int) -> bool:
        """
        Determines if two rings are connected by a bridge based on the list of ring overlaps.

        This method checks if two rings, identified by their IDs, are connected through a bridge by 
        examining the list of ring overlaps. It returns True if a bridge connection is found between the 
        specified rings, and False otherwise.

        Parameters
        :param ring_overlaps: List['RingOverlap']
            A list of RingOverlap objects representing overlaps between rings in the molecule. Each 
            RingOverlap object contains information about the rings involved in the overlap and whether it 
            constitutes a bridge.
        :param ring_id_1: int
            The identifier of the first ring to check for a bridge connection.
        :param ring_id_2: int
            The identifier of the second ring to check for a bridge connection.

        Returns bool:
            True if the specified rings are connected by a bridge, indicating a direct connection that forms 
            part of a bridged ring system, and False otherwise.
        """
        for ring_overlap in ring_overlaps:
            if ring_id_1 == ring_overlap.ring_id_1 and ring_id_2 == ring_overlap.ring_id_2:
                return ring_overlap.is_bridge()
            if ring_id_2 == ring_overlap.ring_id_1 and ring_id_1 == ring_overlap.ring_id_2:
                return ring_overlap.is_bridge()
        return False

    
    def create_bridged_ring(self, involved_ring_ids: Set[int]) -> None:
        """
        Creates a unified representation for a bridged ring system by merging the specified rings into a 
        single RingProperties object.

        This method processes a set of ring IDs that are part of a bridged ring system, creating a new 
        RingProperties object that represents the interconnected rings as a single entity. It involves 
        identifying all atoms and neighbours involved in the bridged system, determining their roles (e.g., 
        bridge atoms), and updating the molecular structure to reflect this unified representation.

        Parameters
        : param involved_ring_ids: Set[int]
            A set of ring IDs that are part of a bridged ring system to be merged into a single 
            RingProperties object.

        Notes
        - Initializes sets for atoms and neighbours involved in the bridged ring system.
        - Iterates through each ring ID in the provided set, marking each as part of a subring of the ridged 
        system and collecting all member atoms and their neighbouring rings.
        - Identifies atoms that are part of the bridged system and classifies them based on their nvolvement 
        in the ring system, distinguishing between those that are bridge atoms and those that are part of he 
        bridged ring itself.
        - Creates a new RingProperties object for the bridged ring, adding it to the class's list of rings 
        and updating its attributes to reflect its bridged nature and interconnectedness.
        - Updates the molecular structure to incorporate the new bridged ring, including updating atom 
        memberships and removing overlaps between the original rings that are now part of the bridged 
        system.
        - This process is crucial for accurately representing complex cyclic structures within the molecule, 
        where rings are interconnected in a way that they share atoms, forming a bridged system. It ensures 
        that the molecular graph accurately reflects the topology of such systems, which is essential for
        the correct calculation of atom positions and the overall layout in 2D space.
        """
        atoms: Set['AtomProperties'] = set()
        neighbours: Set[int] = set()
        for ring_id in involved_ring_ids:
            ring: 'RingProperties' = self.id_to_ring[ring_id]
            for atom in ring.members:
                atoms.add(atom)
            for neighbour_id in ring.neighbouring_rings:
                neighbours.add(neighbour_id)
        leftovers: Set['AtomProperties'] = set()
        ring_members: Set['AtomProperties'] = set()
        for atom in atoms:
            atom_rings_members_id: Set[int] = {ring.id for ring in atom.rings}
            intersect = involved_ring_ids.intersection(atom_rings_members_id)
            if len(atom.rings) == 1 or len(intersect) == 1:
                ring_members.add(atom)
            else:
                leftovers.add(atom)
        for atom in leftovers:
            is_on_ring = False
            for bond, n, m in self.get_bonds_of_atom(atom): 
                atom1, atom2 = self.atoms[n], self.atoms[m] 
                bond_associated_rings = min(len(atom1.rings), len(atom2.rings))
                if bond_associated_rings == 1:
                    is_on_ring = True
            if is_on_ring:
                atom.is_bridge_atom = True
                ring_members.add(atom)
            else:
                atom.is_bridge = True
                ring_members.add(atom)
        bridged_ring = RingProperties(list(ring_members))
        self.add_ring(bridged_ring)
        bridged_ring.bridged = True
        bridged_ring.neighbouring_rings = list(neighbours)
        for ring_id in involved_ring_ids:
            ring = self.id_to_ring[ring_id]
            bridged_ring.subrings.append(ring.copy())
        for atom in ring_members:
            atom.bridged_ring = bridged_ring.id
            for ring_id in involved_ring_ids:
                if self.id_to_ring[ring_id] in atom.rings:
                    atom.rings.remove(self.id_to_ring[ring_id])
            atom.rings.append(bridged_ring)
        involved_ring_ids: List[int] = list(involved_ring_ids)
        for i, ring_id_1 in enumerate(involved_ring_ids):
            for ring_id_2 in involved_ring_ids[i + 1:]:
                self.remove_ring_overlaps_between(ring_id_1, ring_id_2)
        for neighbour_id in neighbours:
            ring_overlaps: List['RingOverlap'] = self.get_ring_overlaps(
                neighbour_id, involved_ring_ids)
            for ring_overlap in ring_overlaps:

                ring_overlap.update_other(bridged_ring.id, neighbour_id)
            neighbour = self.id_to_ring[neighbour_id]
            neighbour.neighbouring_rings.append(bridged_ring.id)


    
    def remove_ring_overlaps_between(self, ring_id_1: int, ring_id_2: int) -> None:
        """
        Removes ring overlaps between two specified rings from the list of ring overlaps.

        This method identifies and removes any ring overlaps between two rings, specified by their IDs, from 
        the class's list of ring overlaps. It ensures that once rings are merged or otherwise processed in a 
        way that eliminates their overlap, the record of their previous overlap is removed to maintain an 
        accurate representation of the molecular structure.

        Parameters
        :param ring_id_1: int
            The identifier of the first ring for which overlaps should be removed.
        :param ring_id_2: int
            The identifier of the second ring for which overlaps should be removed.
        """
        to_remove = []
        for ring_overlap in self.ring_overlaps:
            is_matching_overlap: bool =  \
                (ring_overlap.ring_id_1 == ring_id_1 and ring_overlap.ring_id_2 == ring_id_2) or \
                (ring_overlap.ring_id_2 == ring_id_1 and ring_overlap.ring_id_1 == ring_id_2) 
            if is_matching_overlap:
                to_remove.append(ring_overlap)
                
        for ring_overlap in to_remove:
            self.ring_overlaps.remove(ring_overlap)


    
    def get_ring_overlaps(self, ring_id: int, ring_ids: List[int]) -> List['RingOverlap']:
        """
        Retrieves a list of ring overlaps involving a specified ring and a list of other ring IDs.

        Parameters
        :param ring_id: int
            The identifier of the ring for which overlaps with other rings are to be found.
        :param ring_ids: List[int]
            A list of ring identifiers to check for overlaps with the specified ring.

        Returns List['RingOverlap']
            A list of RingOverlap objects representing the overlaps between the specified ring and any of 
            the rings identified by the IDs in the ring_ids list. Each RingOverlap object contains 
            information about the rings involved in the overlap and the nature of their intersection.
        """
        ring_overlaps: List['RingOverlap'] = []
        for ring_overlap in self.ring_overlaps:
            for ring_id_2 in ring_ids:
                is_matching_overlap: bool = \
                    (ring_overlap.ring_id_1 == ring_id and ring_overlap.ring_id_2 == ring_id_2) or \
                    (ring_overlap.ring_id_2 == ring_id and ring_overlap.ring_id_1 == ring_id_2)
                if is_matching_overlap:
                    ring_overlaps.append(ring_overlap)
        return ring_overlaps
   


    def remove_ring(self, ring: 'RingProperties') -> None:
        """
        Removes a specified ring from the list of rings and updates the list of ring overlaps accordingly.

        Parameters
        :param ring: RingProperties
            The RingProperties object to be removed from the list of rings.
        """
        self.rings.remove(ring)
        overlaps_to_remove = []
        for ring_overlap in self.ring_overlaps:
            if ring_overlap.ring_id_1 == ring.id or ring_overlap.ring_id_2 == ring.id:
                overlaps_to_remove.append(ring_overlap)
        for ring_overlap in overlaps_to_remove:
            self.ring_overlaps.remove(ring_overlap)
        for neighbouring_ring in self.rings:
            if ring.id in neighbouring_ring.neighbouring_rings:
                neighbouring_ring.neighbouring_rings.remove(ring.id)

    def find_bridged_systems(self, rings: List['RingProperties'],
            ring_overlaps: 'RingOverlap') -> List:
        """
        Identifies bridged ring systems within the molecule based on the provided rings and their overlaps.

        Parameters
        :param rings : List['RingProperties']
            A list of RingProperties objects representing the rings within the molecule to be analyzed.
        :param ring_overlaps : List['RingOverlap']
            A list of RingOverlap objects representing overlaps between rings, which is used to determine 
            the interconnectedness of the rings.

        Returns List[List[int]]
            A list of ring groups, where each group is represented as a list of ring IDs. Each group is 
            identified as a bridged system based on the criteria that the number of overlaps is at least as 
            great as the number of rings in the group, indicating a high likelihood of forming a bridged 
            ring system.
        """
        bridged_systems: List = []
        ring_groups = self.get_ring_groups(rings, ring_overlaps)
        for ring_group in ring_groups:
            ring_nr: int = len(ring_group)
            overlap_nr: int = self.get_group_overlap_nr(
                ring_group, ring_overlaps)
            if overlap_nr >= ring_nr:
                bridged_systems.append(ring_group)
        return bridged_systems


    def get_ring_groups(self, rings: List['RingProperties'],
            ring_overlaps: List['RingOverlap']) -> List:
        """
        Organizes rings into groups based on their overlaps, identifying interconnected ring systems within 
        the molecule.

        Parameters
        :param rings: List['RingProperties']
            A list of RingProperties objects representing the rings within the molecule to be analyzed.
        :param ring_overlaps: List['RingOverlap']
            A list of RingOverlap objects representing overlaps between rings, which is used to determine 
            the interconnectedness of the rings.

        Returns List[List[int]]
            A list of ring groups, where each group is represented as a list of ring IDs. Rings within a 
            group are interconnected, either directly or through a series of overlaps, indicating potential 
            bridged or fused ring systems.

        Notes
        - Initializes a list of ring groups, starting with each ring as a separate group.
        - Iteratively merges groups that have overlaps, indicating a structural relationship between rings, 
        until no more merges are possible. This is determined by comparing the number of groups before and 
        after attempting merges.
        - Uses a helper method, ring_groups_have_overlap, to identify if two groups share an overlap, 
        suggesting they should be merged into a single group.
        - Merging is done by creating a union of the two groups and removing the original groups from the 
        list, then adding the merged group. This process simplifies the representation of the molecule's 
        ring structure by consolidating interconnected rings.
        - The merging process continues until the number of groups stabilizes, indicating that all 
        interconnected rings have been grouped together.
        - This method is crucial for simplifying the analysis of molecular structures with complex cyclic 
        components, as it reduces the complexity of the ring structure by grouping interconnected rings. 
        This simplification aids in the identification of bridged and fused ring systems, which are 
        important for accurate layout calculations and visualization.
        - By organizing rings into groups, it provides a basis for further analysis, such as identifying 
        bridged ring systems or resolving the layout of rings in a way that reflects their 
        interconnectedness, which is essential for the accurate representation of molecular topology in 2D 
        space.
        - The final list of ring groups represents a simplified view of the molecule's cyclic structure, 
        where each group may correspond to a bridged, fused, or independent ring system, depending on the 
        overlaps between rings.
        """
        ring_groups = []
        for ring in rings:
            ring_groups.append([ring.id])
        current_ring_nr = 0
        previous_ring_nr = -1
        while current_ring_nr != previous_ring_nr:
            previous_ring_nr = current_ring_nr
            indices = None
            new_group = None
            for i, ring_group_1 in enumerate(ring_groups):
                ring_group_1_found = False
                for j, ring_group_2 in enumerate(ring_groups):
                    if i != j:
                        if self.ring_groups_have_overlap(
                                ring_group_1, ring_group_2, ring_overlaps):
                            indices = [i, j]
                            new_group = list(set(ring_group_1 + ring_group_2))
                            ring_group_1_found = True
                            break
                if ring_group_1_found:
                    break

            if new_group:
                indices.sort(reverse=True)
                for index in indices:
                    ring_groups.pop(index)
                ring_groups.append(new_group)

            current_ring_nr = len(ring_groups)
        return ring_groups

    def ring_groups_have_overlap(self, group_1: List[int], group_2: List[int],
                                 ring_overlaps: List['RingOverlap']) -> bool:
        """
        Determines if two ring groups have an overlap based on the list of ring overlaps.

        Parameters
        :param group_1: List[int]
            The first group of ring IDs to check for overlaps.
        :param group_2: List[int]
            The second group of ring IDs to check for overlaps.
        :param ring_overlaps: List['RingOverlap']
            A list of RingOverlap objects representing overlaps between rings in the molecule. 
            Each RingOverlap object contains information about the rings involved in the overlap.

        Returns bool
            True if an overlap is found between any rings from the two groups, indicating a structural 
            relationship, and False otherwise.
        """
        return any(ring_1 in self.find_neighbouring_rings(ring_overlaps, ring_2)
            for ring_1 in group_1 for ring_2 in group_2)

    @staticmethod
    def get_group_overlap_nr(ring_group,
                             ring_overlaps: List['RingOverlap']) -> int:
        """
        Calculates the number of overlaps within a group of rings based on a list of ring overlaps.

        Parameters:
        :param ring_group List[int]
            A list of ring identifiers (IDs) representing a group of rings to check for overlaps among.
        :param ring_overlaps List['RingOverlap']: 
            A list of `RingOverlap` objects, where each object represents an overlap between two rings,
            identified by their IDs (`ring_id_1` and `ring_id_2`).

        Returns int The total number of overlaps found within the `ring_group`, where an overlap is 
        counted if both rings involved are members of the group.
        """
        ring_group_set = set(ring_group)
        return sum(overlap.ring_id_1 in ring_group_set and overlap.ring_id_2 in ring_group_set
            for overlap in ring_overlaps)

    def get_bonds_of_atom(self, atom: 'AtomProperties') -> Set['Bond']:
        """
        Retrieves all bonds associated with a specified atom within the molecular graph.

        This method searches the molecular graph for bonds that involve a given atom, returning 
        a list of all bonds connected to it. Each bond is represented by a BondProperties 
        object, which encapsulates details about the bond type and the atoms it connects. By 
        iterating through the collection of all bonds in the graph and checking if the specified 
        atom is involved in each bond, the method accurately identifies all connections of the 
        atom, regardless of the atom's role (whether as the starting or ending atom of the 
        bond). 

        Parameters:
        :param atom AtomProperties:
            The atom whose bonds are to be retrieved. This atom is identified by its unique 
            identifier within the molecular graph.

        Returns List[BondProperties]:
            A list of BondProperties objects representing all bonds connected to the specified 
            atom. Each entry in the list corresponds to a distinct bond involving the atom, 
            providing comprehensive information about the atom's connectivity within the 
            molecular structure.
        """
        return set((bond, n, m) for n, m, bond in self.mc.bonds() if atom.id in (n, m))

    def add_ring(self, ring: 'RingProperties') -> None:
        """
        Adds a new ring to the class's collection and updates the internal tracking of ring 
        identifiers.

        Parameters:
        :param ring 'RingProperties': 
            The `RingProperties` object representing the ring to be added. This 
            object encapsulates the properties and characteristics of the ring, such as its member atoms, size, and type (e.g., aromatic, bridged).
        """
        ring.id = self.ring_id_tracker
        self.rings.append(ring)
        self.id_to_ring[ring.id] = ring
        self.ring_id_tracker += 1

    def initial_approximation(self) -> None:
        """
        Determines the initial atom from which to start the layout calculation process for a molecular 
        graph in 2D space.

        This method iterates through the molecular graph to find an appropriate starting atom based on 
        several criteria:
        1. Prefers an atom that is part of a bridged ring system, indicating complex cyclic structures 
        that require careful handling.
        2. If no such atom is found, it looks for an atom that belongs to a bridged ring, which suggests 
        a connection between rings that might need special attention during layout.
        3. If still no suitable atom is found, and if there are rings defined, it selects the first 
        member of the first ring in the class's ring list.
        4. If no rings are defined or none of the above conditions are met, it selects a terminal atom, 
        which is an atom with no more than one bond, simplifying the starting conditions.
        5. As a last resort, if no terminal atom is found, it defaults to the first atom in the graph.

        After selecting the starting atom, it initiates the bond creation process by calling 
        `create_next_bond` with the selected atom, setting the stage for further layout calculations.
        """
        start_atom = None
        for atom in self.graph:
            if atom.bridged_ring is not None:
                start_atom = atom
                break

        if start_atom is None:
            for ring in self.rings:
                if ring.bridged:
                    start_atom = ring.members[0]

        if start_atom is None:
            if len(self.rings) > 0:
                start_ring: 'RingProperties' = self.id_to_ring[0]
                start_atom = start_ring.members[0]

        if start_atom is None:
            for atom in self.graph:
                if atom.is_terminal():
                    start_atom = atom
                    break

        if start_atom is None:
            start_atom = self.graph[0]
        self.create_next_bond(start_atom, None, 0.0)

    def create_next_bond(self, atom: 'AtomProperties',
            previous_atom: Optional['AtomProperties'] = None,
            angle: float = 0.0,
            previous_branch_shortest: bool = False) -> None:
        """
        Creates the next bond for an atom in the molecular structure, updating its position 
        based on the previous atom and angle.

        Parameters:
        :param atom: AtomProperties: 
            The atom for which the next bond is being created.
        :param previous_atom: Optional[AtomProperties]:
            The previous atom connected to the current atom. If None, it is assumed that the 
            current atom is the first in the molecular structure.
        :param angle: float:
            The angle between the previous atom and the current atom in radians. Default is 0.0.
        :param previous_branch_shortest: bool:
            A flag indicating if the previous branch is the shortest. Default is False.

        Logic:
        1. If the atom is already positioned, the method ends without changes.
        2. If there is no previous atom, a special method for the first atom is used.
        3. If the previous atom is connected to one or more rings, a method for calculating the 
        atom's position in ring structures is used.
        4. Otherwise, if the previous atom is not connected to rings, a method for calculating 
        the atom's position without considering rings is used.
        5. If the atom has connected rings, a method for atoms in ring structures is applied.
        6. Depending on the number of neighbors the atom has (from 1 to 4), the corresponding 
        method is chosen to calculate its position, considering various neighbor configurations.
        """
        if atom.positioned:
            return
        if previous_atom is None:
            self.calculate_first_atom(atom)
        elif len(previous_atom.rings) > 0:
            self.calculate_rings(previous_atom, atom)
        else:
            self.calculate_NOT_first_atom(atom, previous_atom, angle)

        if len(atom.rings) > 0:
            self.calculate_some_rings(atom)
        else:
            neighbours: List['AtomProperties'] = atom.neighbours.copy()
            if previous_atom and previous_atom in neighbours:
                neighbours.remove(previous_atom)
            previous_angle: float = atom.get_angle()
            if len(neighbours) == 1:
                self.calculate_1_neighbours(neighbours, atom, previous_atom, previous_angle,
                    previous_branch_shortest)
            elif len(neighbours) == 2:
                self.calculate_2_neighbours(atom, neighbours, previous_atom, previous_angle)
            elif len(neighbours) == 3:
                self.calculate_3_neighbours(atom, neighbours, previous_atom, previous_angle)
            elif len(neighbours) == 4:
                self.calculate_4_neighbours(atom, neighbours, previous_angle)

    def calculate_first_atom(self, atom: 'AtomProperties') -> None:
        """
        Calculates the initial position for the first atom in a molecule.

        This method sets the initial position for the first atom in the molecular structure. It 
        assigns a default position based on the class's bond length and rotates it to a standard 
        orientation. The atom is marked as positioned if it is not part of a bridged ring 
        system, ensuring it's ready for further calculations in the molecular layout process.

        Parameters:
        :param atom AtomProperties:
            The first atom in the molecule to calculate the initial position for.
        """
        dummy: Vector = Vector(self.bond_length, 0)
        dummy.rotate(math.radians(-60.0))
        atom.previous_position = dummy
        atom.previous_atom = None
        atom.set_position(Vector(self.bond_length, 0))
        atom.angle = math.radians(-60.0)
        if atom.bridged_ring is None:
            atom.positioned = True

    def calculate_NOT_first_atom(self, atom: 'AtomProperties',
            previous_atom: 'AtomProperties', angle: float) -> None:
        """
        Calculates the position for an atom that is not the first in the molecule, based on its 
        previous atom and a given angle.

        Parameters:
        :param atom AtomProperties:
            The atom for which the position is being calculated.
        :param previous_atom AtomProperties:
            The atom preceding the current atom in the molecular structure, used as a reference 
            for positioning.
        :param angle float:
            The angle in radians by which the position vector should be rotated to align the 
            atom correctly relative to the previous atom.
        """
        position: Vector = Vector(self.bond_length, 0)
        position.rotate(angle)
        position += previous_atom.position
        atom.set_position(position)
        atom.set_previous_position(previous_atom)
        atom.positioned = True

    def calculate_1_neighbours(self, neighbours: List[int],
            atom: 'AtomProperties', previous_atom: 'AtomProperties',
            previous_angle: float, previous_branch_shortest: bool) -> None:
        """
        Calculates the position for an atom with exactly one neighbor in the molecular 
        structure, considering various bonding scenarios and configurations.

        This method is designed to handle the placement of an atom that has only one neighbor 
        within the molecular structure, taking into account the type of bonds it forms with its 
        previous atom and the presence of any rings. It adjusts the atom's position based on the 
        bond type (single, double, or triple) and the configuration of the molecule, including 
        handling cis and trans isomerism in specific scenarios. The method also considers the 
        angle of the previous bond and the shortest branch condition to correctly orient the 
        atom in space.

        Parameters:
        :param neighbours List[int]:
            A list of atom indices representing the neighbors of the current atom. Since the 
            atom has only one neighbor, this list should contain a single element.
        :param atom AtomProperties:
            The atom for which the position is being calculated.
        :param previous_atom AtomProperties:
            The atom preceding the current atom in the molecular structure, used as a reference 
            for positioning.
        :param previous_angle float:
            The angle in radians between the previous atom and the current atom.
        :param previous_branch_shortest bool:
            Indicates if the previous branch is the shortest, affecting the orientation of the next bond.
        """
        next_atom: 'AtomProperties' = neighbours[0]
        current_bond: Optional['Bond'] = self.bond_lookup(atom, next_atom)
        previous_bond: Optional['Bond'] = None
        if previous_atom:
            previous_bond = self.bond_lookup(previous_atom, atom)
        if current_bond.order == 3 or (previous_bond and previous_bond.order == 3) \
            or (current_bond.order == 2 and previous_bond
            and previous_bond.order == 2 and previous_atom and \
                len(previous_atom.rings) == 0 and len(atom.neighbours) == 2):
            if current_bond.order == 2 and previous_bond.order == 2:
                atom.draw_explicit = True
            if current_bond.order == 3:
                atom.draw_explicit = True
                next_atom.draw_explicit = True
            if current_bond.order == 2 or current_bond.order == 3 or (
                    previous_atom and previous_bond.order == 3):
                next_atom.angle = math.radians(0)
            angle_ = previous_angle + next_atom.angle
            self.create_next_bond(next_atom, atom, angle_)
        elif previous_atom and len(previous_atom.rings) > 0:
            proposed_angle_1: float = math.radians(60.0)
            proposed_angle_2: float = proposed_angle_1 * -1

            proposed_vector_1: Vector = Vector(self.bond_length, 0)
            proposed_vector_2: Vector = Vector(self.bond_length, 0)
            proposed_vector_1.rotate(proposed_angle_1 + atom.get_angle())
            proposed_vector_2.rotate(proposed_angle_2 + atom.get_angle())
            proposed_vector_1 += atom.position
            proposed_vector_2 += atom.position
            centre_of_mass: Vector = self.get_current_centre_of_mass()
            
            distance_1: float = proposed_vector_1.get_squared_distance(centre_of_mass)
            distance_2: float = proposed_vector_2.get_squared_distance(centre_of_mass)
            previous_atom.angle = proposed_angle_2 if distance_1 < distance_2 else proposed_angle_1
            angle_: float = previous_angle + previous_atom.angle
            self.create_next_bond(next_atom, atom, angle_)
        else:
            proposed_angle: float = atom.angle

            if previous_atom and len(previous_atom.neighbours) > 3:
                if round(proposed_angle, 2) > 0.00:
                    proposed_angle: float = min([math.radians(60), proposed_angle])
                elif round(proposed_angle, 2) < 0.00:
                    proposed_angle: float = max([-math.radians(60), proposed_angle])
                else:
                    proposed_angle: float = math.radians(60)
            elif proposed_angle in (0, None):
                last_angled_atom: 'AtomProperties' = self.get_last_atom_with_angle(atom)
                proposed_angle: float = last_angled_atom.angle
                if proposed_angle is None:
                    proposed_angle: float = math.radians(60)

            rotatable: bool = True
            if previous_atom:
                bond: 'Bond' = self.bond_lookup(previous_atom, atom)
                if bond.order == 2:
                    rotatable: bool = False
                    previous_previous_atom: 'AtomProperties' = previous_atom.previous_atom
                    if previous_previous_atom:
                        if (configuration := self.get_configuration(previous_atom, atom)) is not None:
                            if configuration == 'cis':
                                proposed_angle = -proposed_angle
            if rotatable:
                next_atom.angle = proposed_angle if previous_branch_shortest else -proposed_angle
            else:
                next_atom.angle = -proposed_angle
            self.create_next_bond(next_atom, atom,
                previous_angle + next_atom.angle)

    def calculate_2_neighbours(self, atom: 'AtomProperties', neighbours: List['AtomProperties'],
            previous_atom: 'AtomProperties', previous_angle: float) -> None:
        """
        Calculates the positions for an atom with exactly two neighbours in the molecular 
        structure, considering cis and trans isomerism and the shortest branch condition.

        This method is responsible for determining the positions of an atom that has exactly two 
        neighbours within the molecular structure. It takes into account the possibility of cis 
        and trans isomerism and adjusts the atom's orientation based on the shortest branch 
        condition to ensure correct spatial arrangement. The method first checks for the 
        presence of a proposed angle for the atom; if none is found, a default angle is 
        assigned. It then handles cis and trans isomerism by adjusting the angles of the atom 
        and its neighbours accordingly. The method also determines whether the previous branch 
        is the shortest by comparing the sizes of subgraphs involving the previous atom and the 
        neighbours, which influences the orientation of the new bonds created. Finally, it 
        creates the next bonds for the atom with its neighbours, incorporating the calculated 
        angles and the shortest branch condition.

        Parameters:
        :param atom AtomProperties:
            The atom for which the positions are being calculated.
        :param neighbours List[AtomProperties]:
            A list of the atom's neighbours, which should contain exactly two elements.
        :param previous_atom AtomProperties:
            The atom preceding the current atom in the molecular structure, used as a reference 
            for positioning.
        :param previous_angle float:
            The angle in radians between the previous atom and the current atom.
        """
        proposed_angle = atom.angle
        if not proposed_angle:
            proposed_angle = math.radians(60)

        self.handle_cis_trans_isomery(atom, neighbours, previous_atom, proposed_angle)
        subgraph_3_size: int = self.get_subgraph_size(previous_atom, {atom}) if previous_atom else 0

        previous_branch_shortest = False
        if  subgraph_3_size < self.get_subgraph_size(neighbours[0], {atom}) and \
            subgraph_3_size < self.get_subgraph_size(neighbours[1], {atom}):
            previous_branch_shortest = True

        self.create_next_bond(neighbours[0], atom,
            previous_angle + neighbours[0].angle,
            previous_branch_shortest)
        self.create_next_bond(neighbours[1], atom,
            previous_angle + neighbours[1].angle,
            previous_branch_shortest)

    def handle_cis_trans_isomery(self, atom: 'AtomProperties', neighbours: List['AtomProperties'],
            previous_atom: 'AtomProperties', proposed_angle: float) -> None:
        """
        Handles the case of cis and trans isomerism for an atom with two neighbours, adjusting 
        their angles based on the isomeric configuration.

        This method addresses the specific scenario of cis and trans isomerism for an atom 
        connected to two neighbours, determining the correct spatial orientation based on the 
        isomeric configuration and the types of bonds involved. It calculates the subgraph sizes 
        for each neighbour relative to the atom to identify the cis and trans positions, 
        adjusting their angles accordingly. The method also considers the bond types between the 
        atom and its neighbours to further refine the orientation in cases where both bonds are 
        single, potentially adjusting angles based on the configuration of the previous atom in 
        the molecular structure.

        Parameters:
        :param atom AtomProperties:
            The central atom for which cis and trans isomerism is being evaluated.
        :param neighbours List[AtomProperties]:
            A list containing exactly two neighbours of the atom, between which cis and trans 
            isomerism is considered.
        :param previous_atom AtomProperties:
            The atom preceding the current atom in the molecular structure, used for additional 
            configuration checks.
        :param proposed_angle: float: The initial proposed angle for orientation, in radians.
        """
        neighbour_1, neighbour_2 = neighbours
        subgraph_1_size: int = self.get_subgraph_size(neighbour_1, {atom})
        subgraph_2_size: int = self.get_subgraph_size(neighbour_2, {atom})
        cis_atom_index: int = 0
        trans_atom_index: int = 1
        neighbour_1_is_C: bool = ((self.mc.atom(neighbour_1.id).atomic_symbol) == 'C')
        neighbour_2_is_C: bool = ((self.mc.atom(neighbour_2.id).atomic_symbol) == 'C')

        if neighbour_2_is_C and not neighbour_1_is_C and subgraph_2_size > 1 and subgraph_1_size < 5:
            cis_atom_index = 1
            trans_atom_index = 0
        elif not neighbour_2_is_C and neighbour_1_is_C and subgraph_1_size > 1 and subgraph_2_size < 5:
            cis_atom_index = 0
            trans_atom_index = 1
        elif subgraph_2_size > subgraph_1_size:
            cis_atom_index = 1
            trans_atom_index = 0

        cis_atom: 'AtomProperties' = neighbours[cis_atom_index]
        trans_atom: 'AtomProperties' = neighbours[trans_atom_index]

        trans_atom.angle = proposed_angle
        cis_atom.angle = -proposed_angle
        cis_bond: 'Bond' = self.bond_lookup(atom, cis_atom)
        trans_bond: 'Bond' = self.bond_lookup(atom, trans_atom)

        if cis_bond.order == 1 and trans_bond.order == 1:
            if previous_atom:
                previous_bond: 'Bond' = self.bond_lookup(atom, previous_atom)
                if previous_bond.order == 2:
                    if previous_atom.previous_atom:
                        atom1, atom2 = previous_atom, atom
                        configuration_cis_atom: Optional[str] = self.get_configuration(atom1, atom2)
                        if configuration_cis_atom == 'trans':
                            trans_atom.angle = -proposed_angle
                            cis_atom.angle = proposed_angle

    def calculate_3_neighbours(self, atom: 'AtomProperties', neighbours: List['AtomProperties'],
            previous_atom: 'AtomProperties', previous_angle: float) -> None:
        """
        Calculates the positions for an atom with exactly three neighbours in the molecular 
        structure, adjusting angles based on subgraph sizes and ring involvement.

        This method is designed to handle the placement of an atom that has exactly three 
        neighbours within the molecular structure. It determines the orientation of these 
        neighbours based on the sizes of their subgraphs relative to the central atom and 
        adjusts their angles accordingly. The method identifies a 'straight' atom, which is 
        considered the primary direction of extension from the central atom, and two side atoms. 
        The orientation of these atoms is adjusted based on whether they are involved in any 
        rings and the overall structure of the molecule, ensuring a correct spatial arrangement 
        that minimizes overlaps and maintains the integrity of the molecular geometry.

        Parameters:
        :param atom: AtomProperties
            The central atom for which the positions of its neighbours are being calculated.
        :param neighbours List[AtomProperties]:
            A list of the atom's neighbours, which should contain exactly three elements.
        :param previous_atom AtomProperties:
             The atom preceding the current atom in the molecular structure, used as a reference 
             for positioning.
        :param previous_angle float:
             The angle in radians between the previous atom and the current atom, influencing 
             the orientation of the neighbours.
        """
        subgraph_1_size = self.get_subgraph_size(neighbours[0], {atom})
        subgraph_2_size = self.get_subgraph_size(neighbours[1], {atom})
        subgraph_3_size = self.get_subgraph_size(neighbours[2], {atom})
        straight_atom: 'AtomProperties' = neighbours[0]
        left_atom: 'AtomProperties' = neighbours[1]
        right_atom: 'AtomProperties' = neighbours[2]
        if subgraph_2_size > subgraph_1_size and subgraph_2_size > subgraph_3_size:
            straight_atom = neighbours[1]
            left_atom = neighbours[0]
            right_atom = neighbours[2]
        elif subgraph_3_size > subgraph_1_size and subgraph_3_size > subgraph_2_size:
            straight_atom = neighbours[2]
            left_atom = neighbours[0]
            right_atom = neighbours[1]
        if previous_atom and len(previous_atom.rings) < 1\
                and len(straight_atom.rings) < 1\
                and len(left_atom.rings) < 1\
                and len(right_atom.rings) < 1\
                and self.get_subgraph_size(left_atom, {atom}) == 1\
                and self.get_subgraph_size(right_atom, {atom}) == 1\
                and self.get_subgraph_size(straight_atom, {atom}) > 1:
            straight_atom.angle = atom.angle * -1
            if atom.angle >= 0:
                left_atom.angle = math.radians(30)
                right_atom.angle = math.radians(90)
            else:
                left_atom.angle = math.radians(-30)
                right_atom.angle = math.radians(-90)
        else:
            straight_atom.angle = math.radians(0)
            left_atom.angle = math.radians(90)
            right_atom.angle = math.radians(-90)
        self.create_next_bond(straight_atom, atom,
            previous_angle + straight_atom.angle)
        self.create_next_bond(left_atom, atom,
            previous_angle + left_atom.angle)
        self.create_next_bond(right_atom, atom,
            previous_angle + right_atom.angle)

    def calculate_4_neighbours(self, atom: 'AtomProperties',
            neighbours: List['AtomProperties'], previous_angle: float) -> None:
        """
        Handles the case when an atom has exactly four neighbours, adjusting their positions and 
        angles for correct spatial arrangement.

        This method is responsible for calculating the positions and angles of an atom that is 
        connected to four neighbours within the molecular structure. It determines the optimal 
        arrangement of these neighbours based on the sizes of their subgraphs relative to the 
        central atom, ensuring a non-overlapping and structurally sound configuration. The 
        method assigns specific angles to each neighbour to maintain a consistent and clear 
        representation of the molecular geometry, especially in complex molecular structures 
        where an atom is central to four other atoms.

        Parameters:
        :param atom AtomProperties:
            The central atom around which the neighbours are positioned.
        :param neighbours List[AtomProperties]:
            A list of the atom's neighbours, which should contain exactly four elements.
        :param previous_angle float:
            The angle in radians between the previous atom and the current atom, used as a 
            reference for positioning the neighbours.
        """
        subgraph_1_size = self.get_subgraph_size(neighbours[0], {atom})
        subgraph_2_size = self.get_subgraph_size(neighbours[1], {atom})
        subgraph_3_size = self.get_subgraph_size(neighbours[2], {atom})
        subgraph_4_size = self.get_subgraph_size(neighbours[3], {atom})
        atom_1: 'AtomProperties' = neighbours[0]
        atom_2: 'AtomProperties' = neighbours[1]
        atom_3: 'AtomProperties' = neighbours[2]
        atom_4: 'AtomProperties' = neighbours[3]
        if subgraph_2_size > subgraph_1_size and subgraph_2_size > subgraph_3_size\
                and subgraph_2_size > subgraph_4_size:
            atom_1 = neighbours[1]
            atom_2 = neighbours[0]
        elif subgraph_3_size > subgraph_1_size and subgraph_3_size > subgraph_2_size\
                and subgraph_3_size > subgraph_4_size:
            atom_1 = neighbours[2]
            atom_2 = neighbours[0]
            atom_3 = neighbours[1]
        elif subgraph_4_size > subgraph_1_size and subgraph_4_size > subgraph_2_size\
                and subgraph_4_size > subgraph_3_size:
            atom_1 = neighbours[3]
            atom_2 = neighbours[0]
            atom_3 = neighbours[1]
            atom_4 = neighbours[2]

        atom_1.angle = math.radians(-36)
        atom_2.angle = math.radians(36)
        atom_3.angle = math.radians(-108)
        atom_4.angle = math.radians(108)
        self.create_next_bond(atom_1, atom, previous_angle + atom_1.angle)
        self.create_next_bond(atom_2, atom, previous_angle + atom_2.angle)
        self.create_next_bond(atom_3, atom, previous_angle + atom_3.angle)
        self.create_next_bond(atom_4, atom, previous_angle + atom_4.angle)
    


    def get_subgraph_size(self, atom: 'AtomProperties', masked_atoms: Set['AtomProperties']) -> int:
        """
        Calculates and returns the size of a subtree rooted at a given atom, excluding bonds adjacent to atoms specified in masked_atoms.

        This method computes the size of a subtree within a molecular graph, starting from a 
        specified atom and excluding any bonds connected to atoms listed in the `masked_atoms` 
        set. It recursively explores the molecular structure, adding each visited atom to the 
        `masked_atoms` set to avoid revisiting, and counts the total number of unique atoms in 
        the subtree. The size of the subtree is determined by the number of atoms it contains, 
        excluding the initial atom itself.

        Parameters:
        :param atom AtomProperties:
            The root atom from which the subtree size is calculated.
        :param masked_atoms Set[AtomProperties]:
            A set of atoms to be excluded from the subtree calculation, typically used to avoid 
            counting atoms that have already been considered in previous calculations or are not 
            relevant to the current analysis.
        """
        masked_atoms.add(atom)
        for neighbour in atom.neighbours:
            if neighbour not in masked_atoms:
                self.get_subgraph_size(neighbour, masked_atoms)
        return len(masked_atoms) - 1
    

    def calculate_rings(self, previous_atom: 'AtomProperties',
                        atom: 'AtomProperties') -> None:
        """
        Calculates the positions for atoms within rings and aromatic systems, treating the 
        calculation as if it's performed from within the ring itself.

        This method is designed to determine the coordinates of atoms that are part of ring 
        structures or aromatic systems within a molecular graph. It operates under the 
        assumption that the calculation is being performed from the perspective of being inside 
        the ring, allowing for accurate positioning of atoms based on their connectivity and the 
        geometry of the ring. The method takes into account whether the previous atom is part of 
        a bridged ring and adjusts the position of the current atom accordingly, ensuring that 
        the ring's integrity and aromaticity are preserved in the molecular representation. It 
        involves identifying a 'joined vertex' if the previous atom is part of multiple rings, 
        adjusting the position based on the relative positions of neighbours, and setting the 
        atom's position to maintain the ring's structure.

        Parameters:
        :param previous_atom AtomProperties:
             The atom preceding the current atom in the ring, used as a reference for 
             calculating the current atom's position.
        :param atom AtomProperties:
             The atom for which the position is being calculated, ensuring it fits correctly 
             within the ring structure.
        """
        neighbours: List['AtomProperties'] = previous_atom.neighbours
        joined_vertex: Optional['AtomProperties'] = None
        position: Vector = Vector(0, 0)
        if previous_atom.bridged_ring is None and len(previous_atom.rings) > 1:
            for neighbour in neighbours:
                if len(set(neighbour.rings).intersection(set(previous_atom.rings))) == len(previous_atom.rings):
                    joined_vertex: 'AtomProperties' = neighbour
                    break

        if not joined_vertex:
            for neighbour in neighbours:
                if neighbour.positioned and self.atoms_are_in_same_ring(neighbour, previous_atom):
                    position += (neighbour.position - previous_atom.position)
            position.invert()
            position.normalise()
            position *= self.bond_length
            position += previous_atom.position
        else:
            position = joined_vertex.position.copy()
            position.rotate_around_vector(math.pi, previous_atom.position)
        atom.set_previous_position(previous_atom)
        atom.set_position(position)
        atom.positioned = True



    def calculate_some_rings(self, atom: 'AtomProperties') -> None:
        """
        Calculates the coordinates for an atom connected to a ring, handling both bridged and 
        non-bridged ring scenarios.

        This method is responsible for determining the coordinates of an atom that is part of a 
        ring structure within a molecule. It distinguishes between atoms connected to bridged 
        rings and those that are part of regular rings, adjusting the calculation accordingly to 
        ensure accurate positioning within the molecular structure. For atoms connected to 
        bridged rings, it retrieves the specific ring properties to handle the complexity of 
        bridged systems, while for atoms in regular rings, it defaults to the first ring in the 
        atom's ring list. The method then calculates the center position of the ring based on 
        the atom's previous position and the ring's geometry, ensuring the atom is correctly 
        placed relative to the ring's center. This involves inverting the vector from the atom's 
        previous position to its current position, normalizing it, and scaling it according to 
        the ring's radius to find the new center. Finally, it creates the ring with the 
        calculated center, ensuring the atom's position is accurately represented within the 
        ring structure.

        Parameters:
        :param atom AtomProperties:
            The atom for which the ring coordinates are being calculated. This atom is assumed to be part of a ring structure, either directly or through a bridged connection.
        """
        next_ring: 'RingProperties' = self.id_to_ring[atom.bridged_ring] \
            if atom.bridged_ring else atom.rings[0]
            
        if not next_ring.positioned:
            next_center: 'Vector' = atom.previous_position - atom.position
            next_center.invert()            
            next_center.normalise()            
            radius: float = Polygon.find_polygon_radius(self.bond_length, len(next_ring.members))
            next_center *= radius
            next_center += atom.position
            self.create_ring(next_ring, next_center, atom)

    

    def create_ring(self, ring: 'RingProperties', center: Optional[Vector] = None,
            start_atom: Optional['AtomProperties'] = None,
            previous_atom: Optional['AtomProperties'] = None) -> None:
        """
        Creates a ring within a molecular structure, considering its geometry and interaction 
        with other rings.

        This method is responsible for creating and positioning atoms in ring structures of a 
        molecule. It takes into account whether the ring is bridged, determines its center, and 
        arranges atoms according to the geometry of the ring, as well as handles interactions 
        with other rings through common vertices. For bridged rings, a special algorithm is used 
        to correctly position atoms relative to the ring center. For non-bridged rings, atoms 
        are positioned based on a given angle and radius calculated from the bond length and 
        number of members in the ring. The method also handles cases where ring atoms intersect 
        with other rings, ensuring correct connections and preventing overlaps.

        Parameters:
        :param ring RingProperties:
            Properties of the ring for which atom coordinates are being created.
        :param center Optional[Vector]:
            Center of the ring. If not specified, the origin (0, 0) is used.
        :param start_atom Optional[AtomProperties]:
            Starting atom for calculating positions of other atoms in the ring.
        :param previous_atom Optional[AtomProperties]:
            Previous atom used to determine the starting angle.
        """
        if ring.positioned:
            return

        if center is None:
            center: Vector = Vector(0, 0)
        ordered_neighbour_ids: List[int] = self.get_ordered_neighbours(ring, self.ring_overlaps)
        starting_angle: float = 0
        if start_atom:
            starting_angle: float = (start_atom.position - center).angle()
        ring_size: int = len(ring.members)
        radius: float = Polygon.find_polygon_radius(
            self.bond_length, ring_size)
        angle: float = Polygon.get_central_angle(ring_size)
        ring.central_angle = angle
        if start_atom not in ring.members:
            if start_atom:
                start_atom.positioned = False
            start_atom = ring.members[0]

        if ring.bridged:
            KKLayout(structure=self, atoms=ring.members, center=center,
                     start_atom=start_atom, bond_length=self.bond_length)
            ring.positioned = True
            self.set_ring_center(ring)
            center = ring.center
            for subring in ring.subrings:
                self.set_ring_center(subring)
        else:
            self.set_member_positions(ring, start_atom, previous_atom,
                center, starting_angle, radius, angle)
        ring.positioned = True
        ring.center = center

        for neighbour_id in ordered_neighbour_ids:
            neighbour: 'RingProperties' = self.id_to_ring[neighbour_id]
            if neighbour.positioned:
                continue
            atoms: Optional[List['AtomProperties']] = self.get_vertices(
                self.ring_overlaps, ring.id, neighbour.id)
            if len(atoms) == 2:
                self.handle_fused_rings(ring, neighbour, atoms, center)
            elif len(atoms) == 1:
                self.handle_spiro_rings(ring, neighbour, atoms[0], center)

        for atom in ring.members:
            for neighbour in atom.neighbours:
                if neighbour.positioned:
                    continue
                atom.connected_to_ring = True
                self.create_next_bond(neighbour, atom, 0.0)


    def handle_fused_rings(self, ring: 'RingProperties',
            neighbour: 'RingProperties', atoms: List['AtomProperties'],
            center: Vector) -> None:
        """
        Handles the processing of fused cyclic systems within molecular structures, such as 
        decalin ('C12CCCCC1CCCC2').

        This method addresses the specific challenges of dealing with fused ring systems in 
        molecular structures, where two rings share common atoms, creating complex cyclic 
        compounds like decalin. It marks both rings involved as fused, calculates the midpoint 
        between two shared atoms, determines the normals at this midpoint, and adjusts these 
        normals based on the apothem of the neighboring ring to find potential centers for the 
        next ring positions. By comparing distances from a given center to these adjusted 
        normals, it selects the most suitable center for creating a new ring configuration that 
        accommodates the fused nature of the system. Depending on the orientation of the shared 
        atoms, it then creates a new ring with the selected center, ensuring the integrity of 
        the molecular structure is maintained during the fusion process.

        Parameters:
        :param ring RingProperties:
            One of the rings involved in the fusion.
        :param neighbour RingProperties:
            The neighboring ring involved in the fusion.
        :param atoms List[AtomProperties]:
            A list of atoms shared between the fused rings, typically two atoms common to both 
            rings.
        :param center Vector:
            The center point around which the fusion is considered, influencing the orientation 
            of the newly formed ring structure.
        """

        atom_1: 'AtomProperties' = atoms[0]
        atom_2: 'AtomProperties' = atoms[1]
        midpoint: Vector = Vector.get_midpoint(atom_1.position, atom_2.position)
        normals: List[Vector] = Vector.get_normals(atom_1.position, atom_2.position)
        normals[0].normalise()
        normals[1].normalise()

        apothem: float = Polygon.get_apothem_from_side_length(
            self.bond_length, len(neighbour.members))
        normals[0] *= apothem
        normals[0] += midpoint
        normals[1] *= apothem
        normals[1] += midpoint
        next_center: Vector = normals[0]

        distance_to_center_1: float = (center - normals[0]).get_squared_length()
        distance_to_center_2: float = (center - normals[1]).get_squared_length()
        
        if distance_to_center_2 > distance_to_center_1:
            next_center = normals[1]
            
        position_1: Vector = atom_1.position - next_center
        position_2: Vector = atom_2.position - next_center
        
        if position_1.get_clockwise_orientation(position_2) == 'clockwise':
            if not neighbour.positioned:
                self.create_ring(neighbour, next_center, atom_1, atom_2)
        elif not neighbour.positioned:
            self.create_ring(neighbour, next_center, atom_2, atom_1)
        
        
    def handle_spiro_rings(self, ring: 'RingProperties', neighbour: 'RingProperties',
            atom: 'AtomProperties', center: Vector) -> None:
        """
        Handles spirocyclic systems within molecular structures, such as 'C1CCCC11CC1'.

        Spirocyclic systems are characterized by two rings that share a single common atom, 
        forming a spiro junction. This method marks both rings as spirocyclic, calculates a new 
        center for the neighboring ring based on the position of the shared atom and the center 
        of the current ring, and adjusts the position of the neighboring ring accordingly to 
        maintain the integrity of the spirocyclic structure.

        Parameters:
        :param ring RingProperties:
            The current ring involved in the spirocyclic system.
        :param neighbour RingProperties:
            The neighboring ring involved in the spirocyclic system.
        :param atom AtomProperties:
            The atom shared by both rings at the spiro junction.
        :param center Vector:
            The center of the current ring, used as a reference for calculating the new center 
            for the neighboring ring.
        """
        next_center: Vector = center - atom.position
        
        next_center.invert()
        next_center.normalise()
        distance_to_center: float = Polygon.find_polygon_radius(self.bond_length, len(neighbour.members))
        next_center *= distance_to_center
        next_center += atom.position
        if not neighbour.positioned:
            self.create_ring(neighbour, next_center, atom)

    
    def set_ring_center(self, ring: 'RingProperties') -> None:
        """
        Calculates and sets the geometric center of a ring within a molecular structure.

        This method computes the center of a ring by averaging the positions of all atoms that 
        are members of the ring. It iterates through each atom in the ring, summing their 
        positions vectorially, and then divides the total by the number of atoms to find the 
        average position, which represents the ring's center. This center is then assigned to 
        the ring's `center` attribute, providing a reference point for further calculations and 
        manipulations involving the ring.

        Parameters:
        :param ring RingProperties:
             The ring for which the center is to be calculated. This object represents a cyclic 
             structure within the molecule, containing a list of atoms that are part of the ring.
        """
        total: Vector = sum(atom.position for atom in ring.members) / len(ring.members)
        ring.center = total


    def atoms_are_in_same_ring(self, atom_1: 'AtomProperties',
            atom_2: 'AtomProperties') -> bool:
        """
        Determines if two atoms are part of the same ring within a molecular structure.

        This method checks if two given atoms are part of the same ring by comparing their ring 
        memberships. It iterates through the rings of the first atom and checks if any of these 
        rings are also present in the list of rings for the second atom. If any ring is common 
        between the two atoms, it indicates that they are part of the same ring structure.

        Parameters:
        :param atom_1 AtomProperties:
            The first atom to check for ring membership.
        :param atom_2 AtomProperties:
            The second atom to check for ring membership.

        Returns bool:
            True if the atoms are in the same ring, False otherwise.
        """
        return any(ring_id_1 == ring_id_2 for ring_id_1 in atom_1.rings
                   for ring_id_2 in atom_2.rings)

    def get_vertices(self, ring_overlaps: List['RingOverlap'], ring_id_1: int,
                     ring_id_2: int) -> Optional[List['AtomProperties']]:
        """
        Searches for and returns atoms that are in the overlap between two rings identified by ring_id_1 and ring_id_2.

        This method looks for atoms that are located in the overlap between two specified rings. If an overlap is found, it returns a list of atoms that are part of this overlap. If no overlap is found, the method concludes without returning any value, implying a default return of None.

        Parameters:
        :param ring_overlaps List[RingOverlap]: 
            A list of RingOverlap objects representing overlaps between rings in the molecule.
        :param ring_id_1 int:
            The identifier of the first ring to check for overlaps.
        :param ring_id_2 int:
            The identifier of the second ring to check for overlaps.

        Returns Optional[List[AtomProperties]:
            A list of AtomProperties objects representing atoms in the overlap between the two 
            rings, or None if no overlap is found.
        """
        for ring_overlap in ring_overlaps:
            if set((ring_overlap.ring_id_1, ring_overlap.ring_id_2)) == set((ring_id_1, ring_id_2)):
                return [atom for atom in ring_overlap.atoms]


    def get_ordered_neighbours(self, ring: 'RingProperties',
            ring_overlaps: List['RingOverlap']) -> List[int]:
        """
        Retrieves an ordered list of neighboring rings based on the number of atoms they share in common with the specified ring.

        This method is designed to obtain an ordered list of neighboring rings (or structures) 
        based on the number of atoms they have in intersection with the given ring. It returns a 
        list of identifiers for the neighboring rings, sorted by the number of shared atoms in 
        descending order, allowing for the identification of the most closely related rings in 
        terms of atomic overlap.

        Parameters:
        :param ring RingProperties: 
            The ring for which neighboring rings are to be identified and ordered.
        :param ring_overlaps List[RingOverlap]:
            A list of RingOverlap objects representing overlaps between rings in the molecule, 
            used to determine the intersection of atoms between rings.

        Returns List[int]:
            A list of identifiers for neighboring rings, ordered by the number of shared atoms 
            in descending order. Rings with more shared atoms are listed first.
        """
        ordered_neighbours_and_atom_nrs = []
        for neighbour_id in ring.neighbouring_rings:
            atoms: Optional[List['AtomProperties']] = self.get_vertices(
                ring_overlaps, ring.id, neighbour_id)
            ordered_neighbours_and_atom_nrs.append((len(atoms), neighbour_id))

        ordered_neighbours_and_atom_nrs = sorted(
            ordered_neighbours_and_atom_nrs, key=lambda x: x[0], reverse=True)
        ordered_neighbour_ids = [x[1] for x in ordered_neighbours_and_atom_nrs]
        return ordered_neighbour_ids
    

    def set_member_positions(self, ring: 'RingProperties', start_atom: 'AtomProperties',
            previous_atom: Optional['AtomProperties'], center: Vector,
            starting_angle: float, radius: float, angle: float) -> None:
        """
        Positions atoms within a ring structure using polar coordinates (center, radius, angle) 
        and incrementally increases the angle between atoms.

        The primary goal of this method is to arrange atoms in a ring structure by utilizing 
        polar coordinates, where each atom's position is determined relative to a central point, 
        at a specified radius, and with a progressively increasing angle to ensure even 
        distribution around the center. This method iterates through the atoms of the ring, 
        setting their positions based on these polar coordinates until all atoms are positioned 
        or a maximum iteration limit is reached.

        Parameters:
        :param ring RingProperties:
            The ring whose member atoms are to be positioned.
        :param start_atom AtomProperties:
            The starting atom for positioning within the ring.
        :param previous_atom Optional[AtomProperties]:
            The atom preceding the current atom in the positioning sequence. Used as a reference 
            for the first atom's position if it hasn't been positioned yet.
        :param center Vector:
            The central point around which atoms are positioned.
        :param starting_angle float:
            The initial angle in radians for the first atom's position relative to the center.
        :param radius float:
            The distance from the center to the atom's position.
        :param angle float:
            The angular increment in radians between consecutive atoms in the ring.
        """
        current_atom = start_atom
        iteration = 0
        while current_atom is not None and iteration < 100:
            previous = current_atom
            if not previous.positioned:
                x = center[0] + math.cos(starting_angle) * radius
                y = center[1] + math.sin(starting_angle) * radius
                previous.set_position(Vector(x, y))
            starting_angle += angle

            if len(ring.subrings) < 3:
                previous.angle = starting_angle
                previous.positioned = True

            current_atom = self.get_next_in_ring(
                ring, current_atom, previous_atom)
            previous_atom = previous

            if current_atom == start_atom:
                current_atom = None
            iteration += 1


    
    def get_next_in_ring(self, ring: 'RingProperties', current_atom: 'AtomProperties', \
                         previous_atom: 'AtomProperties') -> Optional['AtomProperties']:
        """
        Searches for the next atom in the ring, excluding the previous atom.

        This method iterates through the neighbors of the current atom to find the next atom 
        that is a member of the specified ring, excluding the atom that was previously 
        considered. It checks each neighbor to see if it belongs to the ring and is not the same 
        as the previous atom. If a suitable atom is found, it is returned; otherwise, None is 
        returned. This is useful for traversing the atoms in a ring structure while ensuring 
        that the traversal does not immediately return to the previous atom, allowing for a 
        continuous loop through the ring's members.

        Parameters:
        :param ring RingProperties:
            The ring within which to search for the next atom.
        :param current_atom AtomProperties:
            The current atom from which to start the search.
        :param previous_atom AtomProperties:
             The atom to exclude from the search, typically the atom preceding the current atom 
             in the traversal.

        Returns Optional[AtomProperties]:
            The next atom in the ring that is different from the previous atom, or None if no   
            such atom is found.
        """
        neighbours: List[int] = current_atom.neighbours
        for neighbour in neighbours:
            for member in ring.members:
                if neighbour == member:
                    if previous_atom != neighbour:
                        return neighbour

    @staticmethod
    def find_neighbouring_rings(ring_overlaps: List['RingOverlap'],
            ring_id: int) -> List[int]:
        """
        Finds and returns a list of identifiers for rings neighboring the specified ring.

        This method searches through a list of ring overlaps to identify rings that are adjacent 
        to a given ring, identified by its ID. It examines each overlap to determine if the 
        specified ring is involved and collects the identifiers of neighboring rings, excluding 
        the specified ring itself. This is useful for understanding the connectivity and 
        structure of rings within a molecular graph, especially in complex molecules where rings 
        may overlap or share atoms.

        Parameters:
        :param ring_overlaps List[RingOverlap]:
            A list of RingOverlap objects, each representing an overlap between two rings in the 
            molecule. These objects contain information about which rings are involved in each 
            overlap.
        :param ring_id int:
            The identifier of the ring for which neighboring rings are to be found.

        Returns List[int]:
            A list of identifiers for rings that are neighbors to the specified ring, based on 
            the overlaps. Each identifier represents a ring that shares at least one atom with 
            the specified ring, indicating a direct connection or overlap.
        """
        neighbouring_rings = []
        for ring_overlap in ring_overlaps:
            if ring_overlap.ring_id_1 == ring_id:
                neighbouring_rings.append(ring_overlap.ring_id_2)
            elif ring_overlap.ring_id_2 == ring_id:
                neighbouring_rings.append(ring_overlap.ring_id_1)
        return neighbouring_rings
            

    
    def get_current_centre_of_mass(self) -> Vector:
        """
        Calculates and returns the current center of mass of the molecular graph.

        This method computes the center of mass of the molecular graph by summing the positions 
        of all positioned atoms and dividing by the number of positioned atoms. It iterates 
        through each atom in the graph, adding the position of each atom to a total if the atom 
        is positioned, and then divides this total by the count of positioned atoms to find the 
        average position, which represents the center of mass. This center of mass can be used 
        as a reference point for various calculations and manipulations within the molecular 
        structure, such as aligning or repositioning atoms relative to the overall structure.

        Returns Vector:
            A Vector object representing the coordinates of the center of mass of the molecular 
            graph, based on the positions of all positioned atoms.
        """
        if positioned_atoms := [atom.position for atom in self.graph if atom.positioned]:
            total = np.sum(positioned_atoms, axis=0)
            count = len(positioned_atoms)
            return total / count
        else:
            return Vector(0, 0) 

    def get_last_atom_with_angle(self, atom: 'AtomProperties') -> Optional['AtomProperties']:
        """
        Retrieves the last atom in a chain that has a defined angle relative to the initial atom.

        This method traverses backwards from a given atom through its predecessors until it 
        finds an atom with a defined angle or reaches an atom without a predecessor. It starts 
        with the initial atom's immediate predecessor and continues to trace back through the 
        chain of previous atoms until it encounters an atom with a non-zero angle or reaches the 
        beginning of the chain. 

        Parameters:
        :param atom AtomProperties:
            The starting atom from which to begin the search for an atom with a defined angle.

        Returns Optional[AtomProperties]:
            The last atom in the chain that has a defined angle, or None if no such atom is 
            found before reaching the start of the chain.
        """
        parent_atom: Optional['AtomProperties'] = atom.previous_atom
        angle: float = parent_atom.angle
        while parent_atom and not angle:
            parent_atom = parent_atom.previous_atom
            angle = parent_atom.angle
        return parent_atom

    def get_coord(self, order: List[int]) -> List:
        """
        Packs and returns the coordinates of atoms in a two-dimensional array based on the 
        specified order.

        Parameters:
        :param order List[int]: 
            A list of integers representing the order in which atom coordinates should be 
            packed. Each integer corresponds to an atom index in the class's atom dictionary.

        Returns List[List[float]]:
            A two-dimensional list where each inner list contains the x and y coordinates of an 
            atom, following the order specified in the input list. 
        """
        return [self.atoms[ord].position[:] for ord in order]

    def collision_handling(self) -> None:
        """
        Handles the positioning of all atoms and resolves any overlaps within the molecular structure.

        This method is responsible for adjusting the positions of atoms to minimize overlaps and 
        ensure a visually coherent representation of the molecule. It begins by resolving 
        primary overlaps through a preliminary adjustment phase, followed by an iterative 
        process to refine the positions of atoms based on their bonds and connectivity. The 
        process involves calculating the total overlap score, identifying atoms that can be 
        rotated around their bonds, and adjusting their positions to reduce overlaps. It also 
        considers the complexity of the subgraphs connected to each atom, preferring to rotate 
        smaller subgraphs to minimize disruption. Special handling is given to atoms connected 
        by double bonds, which are less flexible. The method iteratively adjusts the positions 
        of atoms based on their overlap scores, sensitivity thresholds, and the ability to 
        rotate around bonds, aiming to find a configuration that minimizes the total overlap 
        score. Additionally, it accounts for atoms within rings and their specific constraints, 
        applying rotations to subtrees of atoms to resolve overlaps while maintaining the 
        integrity of the molecular structure. The process is repeated for a set number of 
        iterations to gradually improve the layout, with an option for fine-tuning overlaps if 
        enabled.

        The collision handling process includes:
        - Resolving primary overlaps to ensure atoms do not occupy the same space.
        - Calculating the total overlap score and identifying atoms that can be rotated to 
        reduce overlaps.
        - Adjusting atom positions based on their connectivity and the presence of double bonds, 
        which restrict rotation.
        - Considering the depth of subgraphs connected to each atom to decide which atom to 
        rotate in cases of overlap.
        - Applying rotations to subtrees to resolve overlaps, with specific logic for atoms 
        connected by single or double bonds.
        - Repeatedly recalculating overlap scores and adjusting positions until a satisfactory 
        layout is achieved or the maximum number of iterations is reached.
        - Optionally, fine-tuning the positions for further refinement if the `finetune` flag is 
        set.
        """
        self.resolve_primary_overlaps()
        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
        for _ in range(self.overlap_resolution_iterations):
            for atom1_index, atom2_index, bond in self.mc.bonds():
                n: 'AtomProperties' = self.atoms[atom1_index]
                m: 'AtomProperties' = self.atoms[atom2_index]
                if self.can_rotate_around_bond(bond, n, m):
                    tree_depth_1: int = self.get_subgraph_size(n, {m})
                    tree_depth_2: int = self.get_subgraph_size(m, {n})
                    atom_1_rotatable: bool = True
                    atom_2_rotatable: bool = True
                    for neighbouring_bond, *coords in self.get_bonds_of_atom(n):
                        if neighbouring_bond.order == 2:
                            atom_1_rotatable = False
                    for neighbouring_bond, *coords in self.get_bonds_of_atom(m):
                        if neighbouring_bond.order == 2:
                            atom_2_rotatable = False
                    if not atom_1_rotatable and not atom_2_rotatable:
                        continue
                    elif atom_1_rotatable and not atom_2_rotatable:
                        atom_2: 'AtomProperties' = n
                        atom_1: 'AtomProperties' = m
                    elif atom_2_rotatable and not atom_1_rotatable:
                        atom_1: 'AtomProperties' = n
                        atom_2: 'AtomProperties' = m
                    else:
                        atom_1: 'AtomProperties' = m
                        atom_2: 'AtomProperties' = n
                        if tree_depth_1 > tree_depth_2:
                            atom_1: 'AtomProperties' = n
                            atom_2: 'AtomProperties' = m
                    subtree_overlap_score, _ = self.get_subtree_overlap_score(
                        atom_2, atom_1, atom_to_scores)
                    if subtree_overlap_score > self.overlap_sensitivity:
                        neighbours_2 = atom_2.neighbours.copy()
                        neighbours_2.remove(atom_1)

                        if len(neighbours_2) == 1:
                            neighbour = neighbours_2[0]
                            angle = neighbour.position.get_rotation_away_from_vector( \
                                atom_1.position, atom_2.position, math.radians(120))
                            self.rotate_subtree(neighbour, atom_2, angle, atom_2.position)
                            new_overlap_score, _, _ = self.get_overlap_score()
                            if new_overlap_score > self.total_overlap_score:
                                self.rotate_subtree(neighbour, atom_2, -angle, atom_2.position)
                            else:
                                self.total_overlap_score = new_overlap_score
                        elif len(neighbours_2) == 2:
                            if atom_2.rings and atom_1.rings:
                                continue
                            neighbour_1: 'AtomProperties' = neighbours_2[0]
                            neighbour_2: 'AtomProperties' = neighbours_2[1]
                            if len(neighbour_1.rings) == 1 and len(neighbour_2.rings) == 1:
                                if neighbour_1.rings[0] != neighbour_2.rings[0]:
                                    continue
                            elif neighbour_1.rings or neighbour_2.rings:
                                continue
                            else:
                                angle_1 = neighbour_1.position.get_rotation_away_from_vector( \
                                    atom_1.position, atom_2.position, math.radians(120))
                                angle_2 = neighbour_2.position.get_rotation_away_from_vector( \
                                    atom_1.position, atom_2.position, math.radians(120))
                                self.rotate_subtree(neighbour_1, atom_2, angle_1, atom_2.position)
                                self.rotate_subtree(neighbour_2, atom_2, angle_2, atom_2.position)
                                new_overlap_score, _, _ = self.get_overlap_score()
                                if new_overlap_score > self.total_overlap_score:
                                    self.rotate_subtree(neighbour_1, atom_2, -angle_1, atom_2.position)
                                    self.rotate_subtree(neighbour_2, atom_2, -angle_2, atom_2.position)
                                else:
                                    self.total_overlap_score = new_overlap_score
                        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
        for _ in range(self.overlap_resolution_iterations):
            self._finetune_overlap_resolution()
            self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
        for _ in range(self.overlap_resolution_iterations):
            self.resolve_secondary_overlaps(sorted_overlap_scores)

    def resolve_primary_overlaps(self) -> None:
        """
        Resolves initial overlaps in the molecular structure, focusing on cases where a ring has 
        two outgoing edges.

        This method addresses the issue of overlaps that occur when a ring within the molecular 
        structure has two edges extending outwards, which can lead to collisions in the layout, 
        especially noticeable in representations of cyclohexane in a quarter-staggered 
        conformation. It identifies atoms that are part of such overlaps by examining each ring 
        and its members, and then resolves these overlaps by adjusting the positions of the 
        involved atoms. The resolution process involves calculating the angle of rotation needed 
        to minimize the overlap and applying this rotation to the subtrees connected to the 
        overlapping atoms. 
        """
        overlaps: List = []
        resolved_atoms: Dict[int, bool] = {atom.id: False for atom in self.graph}

        for ring in self.rings:
            for atom in ring.members:
                if resolved_atoms[atom.id]:
                    continue
                resolved_atoms[atom.id] = True
                non_ring_neighbours: List['AtomProperties'] = self.get_non_ring_neighbours(
                    atom)
                if len(non_ring_neighbours) > 1 or (len(non_ring_neighbours) == 1 and\
                    len(atom.rings) == 2):
                    overlaps.append({'common': atom, 'rings': atom.rings,
                         'vertices': non_ring_neighbours})
        for overlap in overlaps:
            branches_to_adjust: List['AtomProperties'] = overlap['vertices']
            rings: List['RingProperties'] = overlap['rings']
            root: 'AtomProperties' = overlap['common']
            if len(branches_to_adjust) == 2:
                atom_1, atom_2 = branches_to_adjust
                angle = (2 * math.pi - rings[0].get_angle()) / 6.0
                self.rotate_subtree(atom_1, root, angle, root.position)
                self.rotate_subtree(atom_2, root, -angle, root.position)
                total, sorted_scores, atom_to_score = self.get_overlap_score()
                subtree_overlap_atom_1_1, _ = self.get_subtree_overlap_score(
                    atom_1, root, atom_to_score)
                subtree_overlap_atom_2_1, _ = self.get_subtree_overlap_score(
                    atom_2, root, atom_to_score)
                total_score = subtree_overlap_atom_1_1 + subtree_overlap_atom_2_1
                self.rotate_subtree(atom_1, root, -2.0 * angle, root.position)
                self.rotate_subtree(atom_2, root, 2.0 * angle, root.position)
                total, sorted_scores, atom_to_score = self.get_overlap_score()
                subtree_overlap_atom_1_2, _ = self.get_subtree_overlap_score(
                    atom_1, root, atom_to_score)
                subtree_overlap_atom_2_2, _ = self.get_subtree_overlap_score(
                    atom_2, root, atom_to_score)
                total_score_2 = subtree_overlap_atom_1_2 + subtree_overlap_atom_2_2
                if total_score_2 > total_score:
                    self.rotate_subtree(atom_1, root, 2.0 * angle, root.position)
                    self.rotate_subtree(atom_2, root, -2.0 * angle, root.position)

    @staticmethod
    def get_non_ring_neighbours(atom: 'AtomProperties') -> List['AtomProperties']:
        """
        Identifies and returns a list of neighbours of the specified atom that are not part of 
        any ring it belongs to.

        This method examines the neighbours of a given atom and filters out those that are part 
        of the same rings as the atom, focusing on neighbours that are not involved in any ring 
        structure with the atom. It is useful for understanding the connectivity of an atom 
        within the molecular graph, especially in contexts where the atom's interactions outside 
        of cyclic structures are of interest. By comparing the ring memberships of the atom and 
        its neighbours, it identifies neighbours that do not share any rings with the atom and 
        are not considered bridge atoms, providing insight into the atom's connections to other 
        parts of the molecule that are not part of its immediate cyclic environment.

        Parameters:
        :param atom AtomProperties: 
            The atom for which non-ring neighbours are to be identified.

        Returns: List[AtomProperties]: 
            A list of AtomProperties objects representing neighbours of the specified atom that 
            are not part of any ring the atom is a member of, excluding bridge atoms. 
        """
        non_ring_neighbours: List['AtomProperties'] = []
        for neighbour in atom.neighbours:
            nr_overlapping_rings = len(
                set(atom.ring_indexes).intersection(
                    set(neighbour.ring_indexes)))
            if nr_overlapping_rings == 0 and not neighbour.is_bridge:
                non_ring_neighbours.append(neighbour)
        return non_ring_neighbours

    def rotate_subtree(self, root: 'AtomProperties', root_parent: 'AtomProperties',
            angle: float, center: Vector) -> None:
        """
        Rotates a subtree of the molecular structure around a specified center by a given angle.

        This method rotates a subtree within the molecular graph, starting from a root atom, 
        around a specified center point by a given angle. It is used to adjust the positions of 
        atoms and their associated structures, such as anchored rings, to resolve overlaps or to 
        achieve a desired orientation. The rotation is applied to the root atom and all atoms 
        connected to it within the subtree, excluding the root's parent to maintain the 
        integrity of the molecular structure. This is particularly useful in the layout process 
        to minimize overlaps or to align parts of the molecule according to specific 
        requirements. The rotation affects not only the positions of atoms in the subtree but 
        also updates the centers of any anchored rings associated with the atoms, ensuring that 
        the entire subtree is cohesively repositioned.

        Parameters:
        :param root AtomProperties:
            The root atom of the subtree to be rotated. This atom serves as the starting point 
            for the rotation.
        :param root_parent AtomProperties:
            The parent atom of the root, which is excluded from the rotation to maintain the 
            connection to the rest of the molecular structure.
        :param angle float:
            The angle in radians by which the subtree will be rotated around the center.
        :param center Vector:
            The point around which the rotation is performed. This center is typically the 
            position of a pivotal atom or a calculated point that serves as the axis of rotation.
        """
        for atom in self.traverse_substructure(root, {root_parent}):
            atom.position.rotate_around_vector(angle, center) 
            for anchored_ring in atom.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, center) 

    def traverse_substructure(self, atom: 'AtomProperties',
                              visited: Set['AtomProperties']) -> Generator['AtomProperties', None, None]:
        """
        Traverses a substructure of the molecular graph starting from a given atom, yielding 
        atoms in a depth-first manner.

        This method performs a depth-first traversal of the molecular graph, starting from a 
        specified atom, and yields atoms that are reachable from it, excluding those already 
        visited to avoid cycles. It is a generator function that explores the molecular 
        structure by recursively visiting each atom's neighbours, ensuring that each atom is 
        visited only once. 

        Parameters:
        :param atom AtomProperties:
             The starting atom for the traversal. The traversal begins from this atom and 
             explores the connected substructure.
        :param visited Set[AtomProperties]:
             A set of atoms that have already been visited during the traversal. This is used to 
             avoid revisiting atoms and to ensure that each atom is processed only once.

        Yields AtomProperties:
             Atoms in the connected substructure of the starting atom, yielded one at a time. 
             This allows for iterative processing or analysis of the substructure without the 
             need to construct and return a complete list of atoms upfront, which can be 
             beneficial for large molecular graphs.
        """
        yield atom
        visited.add(atom)
        for neighbour in atom.neighbours:
            if neighbour not in visited:
                yield from self.traverse_substructure(neighbour, visited)


    def get_overlap_score(self) -> Tuple[float, List[Tuple[float, 'AtomProperties']], Dict[int, float]]:
        """
        Calculates the total overlap score and returns a sorted list of atoms by their overlap 
        scores along with the score dictionary.

        This method computes the total overlap score for the molecular graph represented by the 
        class and returns a tuple containing the total overlap score, a list of atoms sorted by 
        their overlap scores in descending order, and a dictionary mapping atom IDs to their 
        individual overlap scores. The overlap score is a measure of how much atoms overlap with 
        each other, indicating the compactness or congestion within the molecular structure. It 
        is calculated based on the distances between atoms, with closer atoms contributing more 
        significantly to the score. The method iterates through all pairs of atoms, calculates 
        the overlap score for each pair based on their distance, and aggregates these scores to 
        determine the total overlap and individual atom overlap scores. The atoms are then 
        sorted by their scores to identify those with the highest overlap, which can be critical 
        for resolving spatial conflicts in the molecular layout.

        Returns Tuple[float, List[Tuple[float, 'AtomProperties'], Dict[int, float]]: 
        A tuple containing:
            - The total overlap score for the molecular graph, representing the overall 
            compactness or congestion.
            - A list of tuples, each containing an atom's overlap score and the atom itself, 
            sorted by the score in descending order. This list helps identify atoms that are 
            most affected by overlaps and may require adjustment.
            - A dictionary mapping atom IDs to their individual overlap scores, providing 
            detailed insights into the distribution of overlaps across the structure.
        """
        total: float = 0.0
        overlap_scores: Dict[int, float] = {atom.id: 0.0 for atom in self.graph}
        atoms: List['AtomProperties'] = list(self.graph)
        positions = np.array([atom.position for atom in atoms])
        distances_squared = np.sum((positions[:, np.newaxis] - positions[np.newaxis, :]) ** 2, axis=-1)
        bond_length_squared = self.bond_length ** 2

        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                if distances_squared[i, j] < bond_length_squared:
                    distance_squared = distances_squared[i, j]
                    weight = (self.bond_length - np.sqrt(distance_squared)) / self.bond_length
                    total += weight
                    overlap_scores[atoms[i].id] += weight
                    overlap_scores[atoms[j].id] += weight
                
        sorted_overlaps: List[Tuple[float, 'AtomProperties']] = \
            [(overlap_scores[atom.id], atom) for atom in atoms]
        sorted_overlaps.sort(key=lambda x: x[0], reverse=True)
        return total, sorted_overlaps, overlap_scores


    def get_subtree_overlap_score(self, root: 'AtomProperties',
                                  root_parent: 'AtomProperties',
                                  atom_to_score: Dict[int, float]) -> Tuple[float, Vector]:
        """
        Calculates the weighted center and total overlap score for a subtree rooted at a given 
        atom, excluding its parent.

        This method computes the total overlap score and the weighted center position for a 
        subtree within the molecular graph, starting from a specified root atom and excluding 
        its parent. The subtree's overlap score is a measure of how much the atoms within the 
        subtree overlap with others, indicating the compactness or congestion of the subtree's 
        layout. The weighted center is calculated based on the positions of atoms that 
        contribute significantly to the overlap, providing a central point that can be used for 
        adjusting the subtree's position to reduce overlaps. The method iterates through the 
        subtree, accumulating the overlap scores of atoms that exceed a sensitivity threshold 
        and adjusting their positions relative to this score to find a central point of gravity. 
        This central point, along with the total score, can guide the repositioning of the 
        subtree to minimize spatial conflicts within the molecular structure.

        Parameters:
        :param root AtomProperties:
            The root atom of the subtree for which the overlap score and weighted center are 
            calculated. This atom serves as the starting point for the traversal and score 
            calculation.
        :param root_parent AtomProperties:
            The parent atom of the root, which is excluded from the subtree to maintain the 
            integrity of the molecular graph's structure during calculations.
        :param atom_to_score Dict[int, float]:
            A dictionary mapping atom IDs to their individual overlap scores, used to determine 
            the contribution of each atom in the subtree to the total overlap score.

        Returns Tuple[float, Vector]: A tuple containing:
            - The average overlap score for the subtree, calculated as the sum of individual 
            atom scores divided by the number of contributing atoms. This score indicates the 
            subtree's overall compactness or the extent of overlap among its atoms.
            - The weighted center position (Vector) of the subtree, derived from the positions 
            of atoms with significant overlap scores. This center is calculated by summing the 
            positions of contributing atoms, each weighted by its overlap score, and then 
            dividing by the total score to find a central point for potential repositioning.
        """
        score = 0.0
        center = Vector(0, 0)
        count = 0
        for atom in self.traverse_substructure(root, {root_parent}):
            subscore = atom_to_score[atom.id]
            if subscore > self.overlap_sensitivity:
                score += subscore
                count += 1
            position = atom.position.copy()
            position *= subscore
            center += position
        if score:
            center /= score
        if count == 0:
            count = 1
        return score / count, center


    def can_rotate_around_bond(self, bond: 'Bond', atom1: 'AtomProperties', atom2: 'AtomProperties') -> bool:
        """
        Determines whether a bond can be rotated to adjust the molecular structure without 
        breaking its integrity.

        This method evaluates whether a given bond within a molecular structure can be rotated 
        as part of layout adjustments, such as resolving overlaps or optimizing the spatial 
        arrangement of atoms. Rotation around a bond is a common operation in molecular graph 
        manipulation, but it's constrained by several factors to maintain the molecule's 
        structural integrity. The method checks the type of the bond, the number of neighbours 
        each atom involved in the bond has, and their involvement in ring structures. 
        Specifically, it assesses whether the bond is a single bond (allowing for rotation), 
        whether either atom has only one neighbour (which would prevent meaningful rotation), 
        and whether both atoms are part of the same ring (which could disrupt the ring's 
        geometry if rotated).

        Parameters:
        :param bond BondProperties:
            The bond to evaluate for potential rotation. This object contains information about 
            the bond type and the atoms it connects.

        Returns bool:
            True if the bond can be safely rotated, indicating that it is a single bond not 
            involving atoms that are solely connected through this bond or are part of the same 
            ring, thereby allowing for adjustments that maintain the molecule's integrity. False 
            otherwise, indicating constraints that prevent rotation to avoid disrupting the 
            molecular structure.
        """
        is_single_bond: bool = (bond.order == 1)
        has_multiple_neighbours: bool = (len(atom1.neighbours) > 1 and len(atom2.neighbours) > 1)
        are_in_same_ring: bool = (atom1.rings and atom2.rings and 
            len(set(atom1.rings).intersection(set(atom2.rings))) > 0)
        return is_single_bond and has_multiple_neighbours and not are_in_same_ring


    def _finetune_overlap_resolution(self) -> None:
        """
        Fine-tunes the resolution of overlaps between atoms in the molecular structure by 
        iteratively adjusting the positions of atoms to minimize the total overlap score.

        This method is designed to refine the positioning of atoms within the molecular 
        structure to reduce overlaps, focusing on atoms that are too close to each other. It 
        operates by identifying pairs of clashing atoms, determining the shortest path between 
        them, and then attempting to rotate the atoms around their bonds to find a configuration 
        that minimizes the overlap. The process involves several steps:

        1. Identifies pairs of atoms that are clashing, i.e., too close to each other based on a 
        predefined sensitivity threshold.
        2. For each pair of clashing atoms, it finds the shortest path connecting them in the 
        molecular graph.
        3. Along this path, it identifies bonds that are rotatable (excluding double bonds) and 
        calculates a distance metric for each bond based on its position in the path.
        4. Selects the bond with the smallest distance metric as the best candidate for rotation 
        to reduce overlap.
        5. Rotates the subtree of atoms around the best bond in increments, evaluating the 
        overlap score after each rotation to find the optimal rotation angle.
        6. Applies the optimal rotation to minimize the total overlap score.

        The method iteratively adjusts the positions of atoms connected by rotatable bonds to 
        resolve overlaps, with a preference for rotating smaller subtrees to minimize structural 
        disruption. It uses a scoring system to evaluate the effectiveness of each rotation and 
        selects the rotation that results in the lowest overlap score. This process is repeated 
        for all identified clashing atom pairs until the total overlap score is below a 
        sensitivity threshold or no further improvement can be made.
        """
        if self.total_overlap_score > self.overlap_sensitivity:
            clashing_atoms: List[Tuple['AtomProperties', 'AtomProperties']] = self._find_clashing_atoms()
            best_connections: List[Tuple[int, int]] = []  # List to store best connections as tuples of atom indices

            for atom_1, atom_2 in clashing_atoms:
                if self.is_connected(atom_1, atom_2):
                    shortest_path: List['AtomProperties'] = self.find_shortest_path(atom_1, atom_2)  
                    rotatable_connections: List[Tuple[int, int]] = []  # Store rotatable connections as tuples
                    distances: List[float] = []
                    for i in range(len(shortest_path) - 1):
                        atom = shortest_path[i]
                        distance_1: int = i
                        distance_2: int = len(shortest_path) - i - 1  # Уменьшаем на 1
                        average_distance = len(shortest_path) / 2
                        distance_metric = abs(average_distance - distance_1) + abs(average_distance - distance_2)

                        if atom.id in self.mc.int_adjacency.keys() and \
                            any(neighbor.id in self.mc.int_adjacency[atom.id] for neighbor in shortest_path):
                            rotatable_connections.append((atom.id, shortest_path[i + 1].id))  # Обращаемся к следующему атома
                            distances.append(distance_metric)

                    best_connection: Optional[Tuple[int, int]] = None
                    optimal_distance: float = float('inf')
                    for i, distance in enumerate(distances):
                        if distance < optimal_distance:
                            best_connection = rotatable_connections[i]
                            optimal_distance = distance

                    if best_connection is not None:
                        best_connections.append(best_connection)

            best_connections = set(best_connections)
            for best_connection in best_connections:
                if self.total_overlap_score > self.overlap_sensitivity:
                    n, m = best_connection  # Unpack the best connection tuple
                    atom_1, atom_2 = self.atoms[n], self.atoms[m] 
                    subtree_size_1: int = self.get_subgraph_size(atom_1, {atom_2})
                    subtree_size_2: int = self.get_subgraph_size(atom_2, {atom_1})
                    if subtree_size_1 < subtree_size_2:
                        rotating_atom = atom_1
                        parent_atom = atom_2
                    else:
                        rotating_atom = atom_2
                        parent_atom = atom_1
                    overlap_score, _, _ = self.get_overlap_score()
                    scores: List[float] = [overlap_score]

                    for i in range(12):
                        self.rotate_subtree(rotating_atom, parent_atom,
                            math.radians(30), parent_atom.position)
                        new_overlap_score, _, _ = self.get_overlap_score()
                        scores.append(new_overlap_score)
                    assert len(scores) == 13
                    scores = scores[:12].copy()
                    best_i = 0
                    best_score = scores[0]
                    for i, score in enumerate(scores):
                        if score < best_score:
                            best_score = score
                            best_i = i
                    self.total_overlap_score = best_score
                    self.rotate_subtree(rotating_atom, parent_atom, \
                        math.radians(30 * best_i + 1), parent_atom.position)


    
    def _find_clashing_atoms(self) -> List[Tuple['AtomProperties', 'AtomProperties']]:
        """
        Identifies and returns a list of atom pairs that are clashing, i.e., positioned too 
        close to each other based on a distance threshold.

        This method scans through all pairs of atoms in the molecular graph to find those that 
        are closer than a specified distance threshold, indicating a clash or overlap. It 
        calculates the squared distance between each pair of atoms and compares it against a 
        threshold value to determine if they are clashing. The threshold is defined as 80% of 
        the squared bond length, aiming to identify atoms that are significantly closer than 
        they should be, considering the typical bond length in the molecular structure. 

        Returns List[Tuple['AtomProperties', 'AtomProperties']:
            A list of tuples, where each tuple contains two 'AtomProperties' objects 
            representing a pair of atoms that are clashing. Each tuple indicates a pair of atoms 
            that are positioned too close to each other, based on the distance threshold.
        """
        clashing_atoms: List[Tuple['AtomProperties', 'AtomProperties']] = []
        atoms: List['AtomProperties'] = list(self.graph.keys())
        for i, atom_1 in enumerate(atoms):
            for j, atom_2 in enumerate(atoms[i + 1:], start=i + 1):
                if self.bond_lookup(atom_1, atom_2) is None:
                    difference: 'Vector' = atom_1.position - atom_2.position
                    distance: float = difference.get_squared_length()
                    if distance < 0.8 * (self.bond_length**2):
                        clashing_atoms.append((atom_1, atom_2))
        return clashing_atoms


    def is_connected(self, atom_1, atom_2) -> bool:
        """
        Determines if two atoms are connected within the molecular graph, i.e., part of the same 
        molecular structure.

        Parameters:
            atom_1 AtomProperties: The first atom to check for connectivity.
            atom_2 AtomProperties: The second atom to check for connectivity.

        Returns bool: 
            True if both atoms are present in the molecular graph, indicating they are part of 
            the same molecular structure, and False otherwise.
        """
        return atom_1 in self.graph and atom_2 in self.graph


    # def find_shortest_path(self, atom_1: 'AtomProperties', atom_2: 'AtomProperties') -> List[Union['Bond', 'AtomProperties']]:
    #     """
    #     Finds the shortest path between two atoms in the molecular graph, returning either the 
    #     sequence of bonds or atoms along the path.

    #     This method implements a shortest path algorithm to determine the most direct route 
    #     between two specified atoms within the molecular structure. It can return the path as a 
    #     list of either the bonds connecting the atoms or the atoms themselves, depending on the 
    #     `path_type` parameter. The algorithm initializes by setting the distance to all atoms as 
    #     infinite, except for the starting atom, which is set to zero. It then iteratively 
    #     selects the atom with the smallest distance that has not been visited, updates the 
    #     distances to its neighbors, and marks it as visited. This process continues until the 
    #     destination atom is reached or all atoms have been visited. Finally, it constructs the 
    #     path from the destination atom back to the starting atom using the recorded previous 
    #     hops.

    #     Parameters
    #     :param atom_1 AtomProperties:
    #         The starting atom from which to find the shortest path.
    #     :param atom_2 AtomProperties:
    #         The destination atom to which to find the shortest path.
    #     :param path_type str:
    #         Specifies the type of elements to return in the path. 'bond' returns the bonds along the path, 'atom' returns the atoms. Default is 'bond'.

    #     Returns List[Union['BondProperties', 'AtomProperties']:
    #         A list representing the shortest path between `atom_1` and `atom_2`. If `path_type` 
    #         is 'bond', the list contains `BondProperties` objects; if 'atom', it contains 
    #         `AtomProperties` objects.

    #     Raises ValueError:
    #         If `path_type` is neither 'bond' nor 'atom'.
    #      """
    #     distances: Dict['AtomProperties', float] = {}
    #     previous_hop: Dict['AtomProperties', Optional['AtomProperties']] = {}
    #     unvisited: Set['AtomProperties'] = set()
    #     for atom in self.graph:
    #         distances[atom] = float('inf')
    #         previous_hop[atom] = None
    #         unvisited.add(atom)
    #     distances[atom_1] = 0.0
    #     while unvisited:
    #         current_atom: Optional['AtomProperties'] = None
    #         minimum: float = float('inf')
    #         for atom in unvisited:
    #             dist: float = distances[atom]
    #             if dist < minimum:
    #                 current_atom: 'AtomProperties' = atom
    #                 minimum = dist
    #         if current_atom is None or current_atom == atom_2:
    #             break
    #         unvisited.remove(current_atom)

    #         for neighbour in self.graph[current_atom]:
    #             if neighbour in unvisited:
    #                 alternative_distance: float = distances[current_atom] + 1.0
    #                 if alternative_distance < distances[neighbour]:
    #                     distances[neighbour] = alternative_distance
    #                     previous_hop[neighbour] = current_atom

    #     path_atoms: List['AtomProperties'] = []
    #     current_atom: Optional['AtomProperties'] = atom_2
    #     if previous_hop[current_atom] or current_atom == atom_1:
    #         while current_atom:
    #             path_atoms.insert(0, current_atom)
    #             current_atom = previous_hop[current_atom]
                
    #     path: List[Union['Bond', 'AtomProperties']] = []
    #     for i in range(1, len(path_atoms)):
    #         atom_1 = path_atoms[i - 1]
    #         atom_2 = path_atoms[i]
    #         bond = self.bond_lookup(atom_1, atom_2)
    #         path.append(bond)
    #     return path 

    def find_shortest_path(self, atom_1: 'AtomProperties', atom_2: 'AtomProperties') -> List['AtomProperties']:
        distances: Dict['AtomProperties', float] = {}
        previous_hop: Dict['AtomProperties', Optional['AtomProperties']] = {}
        unvisited: Set['AtomProperties'] = set()

        # Инициализация расстояний и предыдущих переходов
        for atom in self.graph:
            distances[atom] = float('inf')
            previous_hop[atom] = None 
            unvisited.add(atom)
        
        distances[atom_1] = 0.0
        while unvisited:
            current_atom: Optional['AtomProperties'] = None
            minimum: float = float('inf')
            # Поиск атома с минимальным расстоянием
            for atom in unvisited:
                dist: float = distances[atom]
                if dist < minimum:##
                    current_atom = atom
                    minimum = dist
            # Если нет текущего атома или достигли целевого атома, выходим из цикла
            if current_atom is None or current_atom == atom_2:
                break
            
            unvisited.remove(current_atom)
            # Обновление расстояний до соседей
            for neighbour in self.graph[current_atom]:
                if neighbour in unvisited:
                    alternative_distance: float = distances[current_atom] + 1.0
                    if alternative_distance < distances[neighbour]:
                        distances[neighbour] = alternative_distance
                        previous_hop[neighbour] = current_atom
        # Сбор пути в виде списка атомов
        path_atoms: List['AtomProperties'] = []
        current_atom: Optional['AtomProperties'] = atom_2
        # Проверка на наличие предыдущих переходов и сбор пути
        if previous_hop[current_atom] is not None or current_atom == atom_1:
            while current_atom is not None:
                path_atoms.insert(0, current_atom)
                current_atom = previous_hop[current_atom]
        return path_atoms  # Возвращаем список атомов вместо связей


    def resolve_secondary_overlaps(self, sorted_scores: \
            List[Tuple[float, 'AtomProperties']]) -> None:
        """
        Resolves secondary overlaps in the molecular structure by adjusting the positions of 
        atoms based on their overlap scores.

        This method addresses secondary overlaps in the molecular layout by iteratively 
        adjusting the positions of atoms that have been identified as overlapping beyond a 
        specified sensitivity threshold. It processes a list of atoms sorted by their overlap 
        scores, focusing on atoms with scores higher than a predefined sensitivity level. For 
        each atom, it determines the appropriate action based on the atom's connectivity and 
        proximity to other atoms to minimize overlaps. The process involves finding the closest 
        atom if the target atom has only one neighbor or is isolated, calculating a new position 
        that reduces overlap, and then rotating the atom to this new position. The rotation is 
        designed to move the atom away from its closest neighbor or a specified position to 
        alleviate the overlap.

        Parameters
        :param sorted_scores List[Tuple[float, AtomProperties]]:
            A list of tuples, where each tuple contains an overlap score and an `AtomProperties` 
            object. The list is sorted by score, with the highest scores (indicating more 
            significant overlaps) first.
            
        The steps involved are as follows:
        1. Iterate through atoms sorted by their overlap scores, focusing on those with scores 
        exceeding the sensitivity threshold.
        2. For atoms with one or no neighbors, find the closest atom in the structure to 
        determine a direction for rotation.
        3. Calculate a new position for the atom based on the closest atom's position or a 
        specified reference point, considering the previous positions of both the atom and its 
        closest neighbor.
        4. Rotate the atom to the new position to reduce overlap, using a predefined angle to 
        ensure minimal disruption to the molecular structure.
        """
        for score, atom in sorted_scores:
            if score > self.overlap_sensitivity:
                if len(atom.neighbours) <= 1:
                    if atom.neighbours:
                        continue
                    closest_atom: 'AtomProperties' = self.get_closest_atom(atom)
                    if len(closest_atom.neighbours) <= 1:
                        closest_position: float = closest_atom.previous_position \
                            if closest_atom.previous_position else atom.neighbours[0].position
                    else:
                        closest_position: float = closest_atom.position \
                            if closest_atom.previous_position else atom.neighbours[0].position
                            
                    atom_previous_position: float = atom.previous_position \
                        if atom.previous_position else atom.neighbours[0].position
                        
                    atom.position = atom.position.rotate_away_from_vector(closest_position, \
                        atom_previous_position, math.radians(20))


    
    def get_closest_atom(self, atom: 'AtomProperties') -> 'AtomProperties':
        """
        Identifies and returns the atom closest to a given atom within the molecular graph.

        This method calculates the distance between a specified atom and all other atoms in the 
        molecular graph to find the closest one. It iterates through each atom in the graph, 
        comparing their distances to the given atom and updates the closest atom found so far 
        based on the smallest squared distance. The squared distance is used for efficiency, 
        avoiding the computational cost of square root calculations. The method is useful for 
        determining the spatial relationship between atoms, which can be crucial for layout 
        adjustments, overlap resolution, or identifying nearby atoms for various analyses.

        Parameters
        :param atom AtomProperties:
            The reference atom from which to find the closest atom in the molecular graph.

        Returns AtomProperties:
            The atom determined to be closest to the specified atom based on the shortest 
            distance. If no closer atom is found (e.g., if the graph only contains the specified 
            atom), the method returns None.
        """
        minimal_distance = float('inf')
        closest_atom: Optional['AtomProperties'] = None
        for atom_2 in self.graph:
            if atom == atom_2:
                continue
            squared_distance: float = atom.position.get_squared_distance(atom_2.position)
            if squared_distance < minimal_distance:
                minimal_distance: float = squared_distance
                closest_atom: 'AtomProperties' = atom_2
        return closest_atom


def calculate2d_coord(order, self) -> List[Vector[float, float]]:
        obj = Calculate2d()
        return obj._calculate2d_coord(order, self)

__all__ = ['calculate2d_coord']