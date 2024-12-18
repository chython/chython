"""
This module defines classes that extend the properties of rings, atoms, and bonds within the 
Kaiton structure, focusing on attributes and methods necessary for coordinate calculations.

The classes contained herein serve as supplements to the existing Kaiton structure, offering 
additional functionalities tailored for computational chemistry applications. They are designed 
to facilitate the calculation of molecular geometries by providing detailed attributes and 
methods specific to rings, atoms, and bonds, thereby enhancing the structure's utility in 
algorithms that require precise spatial information. These classes include RingProperties, 
AtomProperties, BondProperties, and RingOverlap, each tailored to represent different aspects of 
molecular structures with attributes and methods that aid in determining spatial relationships 
and characteristics inherent to chemical compounds.

Classes:
- RingProperties: Represents the properties of rings within a molecule, including identifiers, 
    member atoms, positioning status, geometric center, presence of subrings, and types of rings 
    (e.g., bridged, spiro, fused).
- AtomProperties: Encapsulates atomic properties crucial for molecular geometry calculations, 
    such as atomic symbols, positions, and connectivity.
- BondProperties: Details the computational parameters of chemical bonds, including atom 
    references and bond types.
- RingOverlap: Handles overlaps between rings, identifying shared atoms and determining 
    structural characteristics like bridging.

These classes are integral for algorithms that necessitate a deep understanding of molecular 
topology and geometry, offering a structured approach to manipulating and analyzing chemical 
structures programmatically. They facilitate the representation of complex molecular features 
such as ring systems, atomic configurations, and bond characteristics, making them indispensable 
for cheminformatics and computational chemistry applications.
"""

from typing import List, Optional, Tuple
from .MathHelper import Vector
import math

class RingProperties:
    """
    A class on computing parameters of rings
    """
    def __init__(self: 'RingProperties', ring: List['AtomProperties']) -> None:
        """
        Constructor of the class that complements information about rings in the Kaiton 
        structure, creating new properties or converting existing ones into a more convenient 
        form for use in coordinate calculation algorithms.

        Parameters:
        :param ring List['AtomProperties']: 
            A list of AtomProperties objects forming the current ring.

        Attributes:
        - id Optional[int]: Contains information about the ring identifier, its sequential 
            number.
        - members List['AtomProperties']: A list of references to corresponding AtomProperties 
            objects which are participants of this ring.
        - members_id List[int]: A list of identifiers for AtomProperties objects which are 
            participants of this ring.
        - positioned bool: A boolean value corresponding to the state of the ring calculation, 
            returns True if all atoms of this ring have received their coordinates.
        - center 'Vector': The center of the ring in coordinates.
        - subrings List: A boolean value indicating whether the ring contains subrings.
        - bridged bool: A boolean value indicating whether the ring is a bridge ring.
        - spiro bool: A boolean value indicating whether the ring is a spirocycle.
        - fused bool: A boolean value indicating whether it is part of a condensed cyclic system.
        - subring_of_bridged bool: A boolean value indicating whether the subrings are bridged.
        - central_angle float: The central angle of the ring.
        - neighbouring_rings List[int]: A list of identifiers of neighboring rings to the 
            current one.
        """
        self.id: Optional[int] = None
        self.members: List['AtomProperties'] = ring
        self.members_id: List[int] = [atom.id for atom in self.members]
        self.positioned: bool = False
        self.center: 'Vector' = Vector(0, 0)
        self.subrings: List = []
        self.bridged = False
        self.spiro: bool = False
        self.fused: bool = False
        self.subring_of_bridged = False
        self.central_angle: float = 0.0
        self.neighbouring_rings: List[int] = [] 

        # добавляем в свойства атомов то что они находятся в этом кольце
        for atom in self.members:
            atom.ring_indexes.append(self.id)
            atom.rings.append(self)


    def __eq__(self, other: 'RingProperties') -> bool:
        """
        Compares two RingProperties instances for equality based on their identifiers.

        This method checks if the identifiers of the two RingProperties instances being compared 
        are equal, returning True if they match and False otherwise. It serves as a quick way to 
        determine if two rings refer to the same entity in terms of their unique identifier.

        Parameters:
        :param other RingProperties:
            The instance to compare with the current instance.

        Returns bool: 
            True if the identifiers of the two instances are equal, indicating they represent 
            the same ring. False otherwise.
        """
        return False if other is None else self.id == other.id
    
    
    def __hash__(self) -> int:
        """
        Returns the hash value of the current object, which is the unique identifier of the ring.

        This method is used when the object needs to be inserted into a hash-based collection 
        such as a set or dictionary.
        The hash value is derived from the ring's unique identifier, allowing for efficient 
        storage and retrieval of ring objects
        in collections that rely on hashing.

        Returns int: 
            The unique identifier of the current object of the class, used as the hash value.
        """
        return self.id
    
    
    def get_angle(self) -> float:
        """
        Calculates the exterior angle of the polygon formed by the ring in radians.

        The exterior angle is determined by subtracting the central angle of the ring from π 
        (pi), providing a measure
        of the angle formed outside the ring by extending one of its sides. This method is 
        particularly useful for
        understanding the geometry of the ring within the context of its surrounding environment.

        Returns float: 
            The exterior angle of the polygon formed by the ring in radians.
        """
        return math.pi - self.central_angle
    
    
    def __repr__(self) -> str:
        """
        Provides a human-readable representation of the RingProperties object, primarily 
        intended for debugging purposes.

        The representation includes the ring's identifier followed by the identifiers of the 
        atoms that are part of the ring, separated by hyphens. This format offers a concise yet 
        informative overview of the ring's composition, aiding in the identification and 
        analysis of rings during development and debugging sessions.

        Returns str: 
            A string combining the ring's identifier and the identifiers of the atoms that make 
            up the ring, separated by hyphens.
        """
        members: str = '-'.join(str(member) for member in self.members)
        return f'{self.id} {members}'


    def copy(self) -> 'RingProperties':
        """
        Creates a deep copy of the current RingProperties instance, duplicating all its 
        attributes and relationships.

        This method constructs a new RingProperties object that mirrors the current 
        instance exactly, including the list of member atoms, their identifiers, the 
        ring's position status, geometric center, presence of subrings,
        and various boolean flags indicating the ring's characteristics (e.g., whether it 
        is bridged, spiro, fused).
        Additionally, it copies over the list of neighboring rings and any subrings 
        associated with the ring.

        Returns RingProperties: 
            A new instance of the RingProperties class that is a deep copy of the current 
            instance, complete with all attributes and relationships duplicated.
        """
        new_members: List['AtomProperties'] = []
        for atom in self.members:
            new_members.append(atom.copy())

        new_ring = RingProperties(new_members)
        new_ring.id = self.id
        for ring_id in self.neighbouring_rings:
            new_ring.neighbouring_rings.append(ring_id)

        new_ring.positioned = self.positioned
        for subring in self.subrings:
            new_ring.subrings.append(subring)
        new_ring.bridged = self.bridged
        new_ring.subring_of_bridged = self.subring_of_bridged
        new_ring.spiro = self.spiro
        new_ring.fused = self.fused
        new_ring.central_angle = self.central_angle
        return new_ring
    























class BondProperties:
    """
    A class about computational parameters of links
    """
    def __init__(self: 'BondProperties', atom1: 'AtomProperties', \
                 atom2: 'AtomProperties', bond) -> None:
        """
        Constructor of the class that complements information about bonds in the Kaiton 
        structure, creating new properties or converting existing ones into a more 
        convenient form for use in coordinate calculation algorithms.

        Parameters:
        :param atom1 'AtomProperties': 
            Reference to the object of the class of the first atom forming this bond.
        :param atom2 'AtomProperties': 
            Reference to the object of the class of the second atom forming this bond.
        :param bond: 
            Reference to the original Kaiton bond class.

        Attributes:
        - id Tuple['AtomProperties']: Identifier of the current bond, which is a tuple of 
            atom identifiers between which this bond exists.
        - n int: Identifier of the first atom of this bond.
        - m int: Identifier of the second atom of this bond.
        - atom1 'AtomProperties': Reference to the object of the class of the first atom of 
            this bond.
        - atom2 'AtomProperties': Reference to the object of the class of the second atom of 
            this bond.
        - type str: String that characterizes the type of bond, primary, secondary, or 
            tertiary.

        # center (bool): Placeholder for future expansion.
        # chiral (bool): Placeholder for future expansion.
        # chiral_symbol (Optional[str]): Placeholder for future expansion.

        The constructor initializes the bond properties based on the provided atoms and 
        determines its type (single, double, triple) based on the order of the bond.
        
        """
        self.id: Tuple['AtomProperties'] = (atom1.id, atom2.id)
        self.n: int = atom1.id #atom1 index
        self.m: int = atom2.id #atom2 index

        self.atom1: 'AtomProperties' = atom1
        self.atom2: 'AtomProperties' = atom2

        # self.center: bool = False # рудименты кода
        # self.chiral: bool = False # рудименты кода
        # self.chiral_symbol: Optional[str] = None # # рудименты кода

        self.type = Optional[None]
        if bond.order == 1:
            self.type = 'single'
        elif bond.order == 2:
            self.type = 'double'
        elif bond.order == 3:
            self.type = 'triple'


























class AtomProperties:
    """
    A class about computing parameters of atoms
    """
    def __init__(self: 'AtomProperties', atom_index: int, symbol: str) -> None:
        """
        Initializes an instance of the AtomProperties class with data about an atom.

        Parameters:
        :param atom_index int: 
            The index of the current atom within the molecular structure.
        :param symbol str: 
            Symbol representing the element according to the periodic table.

        Attributes:
        - id int: Unique identifier for the atom.
        - symbol str: String characterizing the name of the element according to the periodic 
        table.
        - ring_indexes List[int]: List of identifiers for rings in which the atom is a 
        participant.
        - rings List['RingProperties']: List of references to RingProperties objects 
        representing the rings in which the atom is involved.
        - is_bridge_atom bool: Flag indicating whether the atom is a bridging atom.
        - is_bridge bool: Flag indicating whether the atom is a bridging atom.
        - bridged_ring Optional['RingProperties']: RingProperties object representing the ring 
        through which a bridge passes.
        - positioned bool: Flag indicating whether the coordinates for the current atom have 
        been calculated.
        - previous_position 'Vector': Coordinates of the preceding atom.
        - position 'Vector': Current coordinates of the atom.
        - angle Optional[float]: Angle between the current and preceding atom.
        - force_positioned bool: Flag indicating whether the atom's position was calculated 
        forcibly.
        - connected_to_ring bool: Flag indicating whether the atom is connected to a ring.
        - draw_explicit bool: Flag indicating whether the atom should be drawn explicitly.
        - neighbours List['AtomProperties']: List of AtomProperties objects representing atoms 
        with which the current atom forms bonds.
        - previous_atom Optional['AtomProperties']: Reference to an AtomProperties object 
        representing the preceding atom in the chain.
        """

        self.id: int = atom_index
        self.symbol: str = symbol
        self.ring_indexes: List[int] = []
        self.rings: List['RingProperties'] = []

        # self.original_rings: List['RingProperties'] = []
        self.anchored_rings: List['RingProperties'] = []
        self.is_bridge_atom: bool = False
        self.is_bridge: bool = False

        self.bridged_ring = None
        self.positioned: bool = False

        self.previous_position: 'Vector' = Vector(0, 0)
        self.position: 'Vector' = Vector(0, 0)
        self.angle: Optional[float] = None
        self.force_positioned: bool = False
        self.connected_to_ring: bool = False
        self.draw_explicit: bool = False
        self.neighbours: List['AtomProperties'] = []
        self.previous_atom: Optional['AtomProperties'] = None

    
    def __eq__(self, other: 'AtomProperties') -> bool:
        """
        Compares two AtomProperties instances for equality based on their identifiers.

        Parameters:
        :param other 'AtomProperties': 
            Another instance of the AtomProperties class.

        Returns bool:
            True if both instances represent atoms with the same identifier, otherwise False.
        """
        return False if other is None else self.id == other.id
    
    
    def set_position(self, vector: 'Vector') -> None:
        """
        Sets the position of the current atom to the specified vector.

        Parameters:
        :param vector 'Vector': 
            An instance of the Vector class, whose coordinates are assigned as the position of the current atom.
        """
        self.position: 'Vector' = vector

    
    def __hash__(self) -> int:
        """
        Returns the hash value of the current object, which is the unique identifier of the atom.

        This method is used when the object needs to be inserted into a hash-based collection 
        such as a set or dictionary.

        Returns int:
            The unique identifier of the current object of the class, used as the hash value.
        """
        return self.id
    
    
    def __repr__(self) -> str:
        """
        Provides a human-readable representation of the AtomProperties object, primarily 
        intended for debugging purposes.

        The representation includes the atomic symbol followed by the atomic index, adjusted by 
        subtracting 1 due to the indexing convention in Chython where numbering starts from 1 
        instead of 0.

        Returns str:
            A string combining the atomic symbol and the adjusted atomic index.
        """
        return f'{self.symbol}_{self.id - 1}'
    
    
    def get_angle(self, reference_vector: Optional['Vector']=None) -> float:
        """
        Calculates the angle between the current atom and either the previous atom or a 
        specified reference vector.

        By default, the angle is calculated between the current atom and the previous atom. 
        However, if a reference_vector is provided, the angle between the current atom and the 
        reference_vector will be calculated instead.

        Parameters:
        :param reference_vector Optional['Vector']: 
            An object of the Vector class representing the coordinates with which the angle will 
            be calculated. If None, the angle between the current atom and the previous atom is 
            calculated. Defaults to None.

        Returns float: 
            The angle between the current atom and either the previous atom or the specified 
            reference vector, depending on the parameter provided.
        """
        vector_1: float = self.position
        vector_2: float = self.previous_position if not reference_vector else reference_vector
        vector = Vector.subtract_vectors(vector_1, vector_2)
        return vector.angle()

    
    def copy(self) -> 'AtomProperties':
        """
        Creates a deep copy of the current AtomProperties instance and returns it as a new 
        object of the same class.

        This method duplicates all attributes of the current atom, including its position, 
        connections, and identifiers, ensuring that modifications to the copy do not affect the 
        original atom object.

        Returns AtomProperties: 
            A new instance of the AtomProperties class with identical properties to the original 
            atom, but as a separate object in memory.
        """
        new_atom = AtomProperties(self.id, self.symbol)
        new_atom.ring_indexes =self.ring_indexes
        new_atom.rings = self.rings
        # new_atom.original_rings = self.original_rings
        new_atom.anchored_rings = self.anchored_rings
        new_atom.is_bridge_atom = self.is_bridge_atom
        new_atom.is_bridge = self.is_bridge
        new_atom.positioned = self.positioned
        new_atom.previous_position = self.previous_position
        new_atom.position = self.position
        new_atom.angle = self.angle
        new_atom.force_positioned = self.force_positioned
        new_atom.connected_to_ring = self.connected_to_ring
        new_atom.draw_explicit = self.draw_explicit
        new_atom.neighbours = self.neighbours
        new_atom.previous_atom = self.previous_atom
        return new_atom
    

    def is_terminal(self) -> bool:
        "Returns boolean whether a given atom is terminal (has no more than one bond)."
        return len(self.neighbours) <= 1

    def set_previous_position(self, previous_atom: 'AtomProperties') -> None:
        "Set previous position atom"
        self.previous_position = previous_atom.position
        self.previous_atom = previous_atom













class RingOverlap:
    """
    Initializes an instance of the RingOverlap class, which represents the overlap between 
    two rings.
    """
    def __init__(self, ring_1: 'RingProperties', ring_2: 'RingProperties') -> None:
        """
        This class is designed to handle situations where two rings share common atoms, 
        indicating an overlap or intersection between them.

        Parameters:
        :param ring_1 RingProperties: 
            An instance of the RingProperties class representing the first ring involved in the 
            overlap.
        :param ring_2 RingProperties: 
            An instance of the RingProperties class representing the second ring involved in the 
            overlap.

        Attributes:
        - id: A unique identifier for the overlap instance. Initially set to None, indicating 
        that the overlap ID may need to be assigned externally.
        - ring_id_1 int: The identifier of the first ring participating in the overlap.
        - ring_id_2 int: The identifier of the second ring participating in the overlap.
        - atoms List['AtomProperties']: A list of AtomProperties instances that are common to 
        both rings, representing the atoms where the overlap occurs.

        The constructor identifies the common atoms between the two rings by intersecting the 
        members of both rings and stores their identifiers for reference.
        """
        self.id = None
        self.ring_id_1: int = ring_1.id
        self.ring_id_2: int = ring_2.id
        self.atoms: List['AtomProperties'] = set(ring_1.members).intersection(set(ring_2.members))


    def __repr__(self) -> str:
        """
        Provides a human-readable representation of the RingOverlap object, primarily intended 
        for debugging purposes.

        Returns a string that includes the overlap identifier and the identifiers of the rings 
        involved in the overlap, offering a convenient way to inspect the state of the object 
        quickly.

        Returns str: 
            A formatted string containing the overlap's unique identifier (`id`) and the 
        identifiers of the two rings (`ring_id_1` and `ring_id_2`) participating in the overlap. 
        This representation aids in quickly identifying the instance's state, especially useful 
        during debugging sessions.
        """
        return f'{self.id=}, {self.ring_id_1=}, {self.ring_id_2=}'


    def is_bridge(self) -> bool:
        """
        Determines whether the overlap represents a bridge between rings based on the number of 
        common atoms and their ring memberships.

        This method checks if the current overlap involves more than two atoms or if any atom 
        within the overlap participates in more than two rings, indicating a bridging structure.

        Returns bool: 
            True if the overlap is considered part of a single bridge ring, otherwise False. 
            Specifically, returns True if either the number of atoms involved in the overlap 
            exceeds two or if any atom in the overlap belongs to more than two rings, suggesting 
            a complex bridging configuration.
        """
        return len(self.atoms) > 2 or any(len(atom.rings) > 2 for atom in self.atoms)
        # ниже старая версия функции
        # if len(self.atoms) > 2:
        #     return True
        # for atom in self.atoms:
        #     if len(atom.rings) > 2:
        #         return True
        # return False


    def involves_ring(self, ring_id: int) -> bool:
        """
        Checks if a ring identified by `ring_id` is involved in the current overlap.

        This method determines whether the specified ring identifier matches either of the two 
        rings participating in the overlap represented by the current instance of RingOverlap.

        Parameters:
        :param ring_id int: 
            The identifier of the ring to check for involvement in the overlap.

        Returns bool: 
            True if the specified ring identifier matches either `ring_id_1` or `ring_id_2`, 
            indicating that the ring is part of the current overlap. False otherwise.
        """
        return self.ring_id_1 == ring_id or self.ring_id_2 == ring_id
    

    def update_other(self, ring_id: int, other_ring_id: int) -> None:
        """
        Updates the current attributes of the class depending on the other ring.

        This method adjusts the internal identifiers of the RingOverlap instance based on the 
        provided ring identifiers, ensuring that the instance accurately reflects the rings 
        involved in the overlap.

        Parameters:
        :param ring_id int: 
            Identifier of the first ring to be considered for updating.
        :param other_ring_id int: 
            Identifier of another ring, used to determine which attribute (ring_id_1 or 
            ring_id_2) needs to be updated.

        If the current instance's ring_id_1 matches other_ring_id, then ring_id is assigned to 
        ring_id_2, and vice versa. This ensures that the RingOverlap instance correctly tracks 
        the two rings involved in the overlap.
        """
        if self.ring_id_1 == other_ring_id:
            self.ring_id_2 = ring_id
        else:
            self.ring_id_1 = ring_id