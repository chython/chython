"""
Defines the KKLayout class, utilized for arranging molecular structures using the Kamada-Kawai 
algorithm.

Description:
The KKLayout class is designed to optimize the layout of molecular structures through the 
application of the Kamada-Kawai algorithm, a graph drawing method that models atoms as physical 
bodies connected by springs, aiming to find a low-energy configuration. This approach allows for 
the visualization of molecular structures in a two-dimensional plane, where atoms are positioned 
according to calculated coordinates to reflect their interconnections and minimize the system's 
total energy. The class is initialized with essential details about the molecular structure, 
enabling the creation of matrices for storing atomic distances, bond stiffness, and interaction 
energies. These matrices support the iterative process of finding an optimal layout that 
minimizes energy, thereby enhancing the stability and clarity of the molecular representation.

Outcome:
Achieves a stable molecular layout where atomic positions are optimized for minimal energy, 
enhancing the interpretability of complex molecular structures. This optimized arrangement not 
only reflects the underlying chemistry but also simplifies the visual analysis of molecular 
architecture, making it invaluable for scientific and educational purposes.
"""
from typing import Dict, List, TYPE_CHECKING, Tuple
from .MathHelper import Vector, Polygon
import math
if TYPE_CHECKING:
    from .Calculate2d import Calculate2d
    from .Properties import *

class KKLayout:
    """
    Class for calculating the optimal arrangement of atoms in a molecular structure
    using the Kamada-Kawai algorithm.
    """
    def __init__(self, structure: 'Calculate2d', atoms: List['AtomProperties'], \
                 center: 'Vector', start_atom: 'AtomProperties', bond_length: float, \
                 threshold: float=0.1, inner_threshold: float=0.1, max_iteration: int=2000,
                 max_inner_iteration: int=50, max_energy: int=1e9):
        """
        Initializes the KKLayout object with the necessary parameters to calculate the optimal
        arrangement of atoms in a molecular structure using the Kamada-Kawai algorithm.

        Parameters:
        :param structure Calculate2d: 
            An instance of the Calculate2d class used for performing 2D space calculations.
        :param atoms List[AtomProperties]: 
            A list containing instances of the AtomProperties class representing the atoms
            in the molecular structure.
        :param center Vector: 
            An instance of the Vector class specifying the geometric center of the molecular 
            structure.
        :param start_atom AtomProperties: 
            An instance of the AtomProperties class designating the starting atom for 
            constructing the layout.
        :param bond_length float: 
            The specified length of bonds between atoms in the molecular structure.
        :param threshold Optional[float]: 
            The energy threshold value used to determine when to stop iterations during the 
            calculation. Defaults to 0.1.
        :param inner_threshold Optional[float]: 
            The inner energy threshold value used to determine when to stop internal iterations
            during the calculation. Defaults to 0.1.
        :param max_iteration Optional[int]: 
            The maximum number of iterations allowed for the calculation process. Defaults to 
            2000.
        :param max_inner_iteration Optional[int]: 
            The maximum number of inner iterations allowed for the calculation process.
            Defaults to 50.
        :param max_energy Optional[int]: 
            The maximum allowable energy level for the system being calculated. Defaults to 1e9.

        Attributes:
        self.structure: Stores the Calculate2d instance passed as a parameter.
        self.atoms: Stores the list of AtomProperties instances representing the atoms in the 
            structure.
        self.center: Stores the Vector instance representing the center of the molecular 
            structure.
        self.start_atom: Stores the AtomProperties instance representing the starting atom for 
            the layout construction.
        self.edge_strength: Stores the bond length between atoms.
        self.threshold: Stores the energy threshold value for stopping iterations.
        self.inner_threshold: Stores the inner energy threshold value for internal iterations.
        self.max_iteration: Stores the maximum number of iterations allowed.
        self.max_inner_iteration: Stores the maximum number of inner iterations allowed.
        self.max_energy: Stores the maximum allowable energy level for the system.

        Additional Attributes:
        self.x_positions, self.y_positions: Dictionaries storing the X and Y coordinates of each 
            atom in the structure.
        self.positioned: A dictionary indicating whether each atom has been positioned.
        self.length_matrix: A matrix storing the lengths of bonds between pairs of atoms.
        self.distance_matrix: A matrix storing the distances between pairs of atoms.
        self.spring_strengths: A matrix storing the stiffness values of links between pairs of 
            atoms.
        self.energy_matrix: A matrix storing the interaction energy between pairs of atoms.
        self.energy_sums_x, self.energy_sums_y: Dictionaries storing the summations of energies 
            along the X and Y axes for each atom.
        """
        self.structure: 'Calculate2d' = structure
        self.atoms: List['AtomProperties'] = atoms
        self.center: 'Vector' = center
        self.start_atom: 'AtomProperties' = start_atom
        self.edge_strength: int = bond_length
        self.threshold: float = threshold
        self.inner_threshold: float = inner_threshold
        self.max_iteration: int = max_iteration
        self.max_inner_iteration: int = max_inner_iteration
        self.max_energy: int = max_energy
        
        self.x_positions: Dict['AtomProperties', float] = {}
        self.y_positions: Dict['AtomProperties', float] = {}
        self.positioned: Dict['AtomProperties', bool] = {}
        self.length_matrix: Dict['AtomProperties', Dict['AtomProperties', float]] = {}
        self.distance_matrix: Dict[int, Dict[int, float]] = {}
        self.spring_strengths: Dict['AtomProperties', Dict['AtomProperties', float]] = {}
        self.energy_matrix: Dict['AtomProperties', Dict['AtomProperties', Optional[float]]] = {}
        self.energy_sums_x: Dict['AtomProperties', Dict['AtomProperties', Optional[float]]] = {}
        self.energy_sums_y: Dict['AtomProperties', Dict['AtomProperties', Optional[float]]] = {}
        
        self.initialise_matrices()
        self.get_kk_layout()
  
  
    def initialise_matrices(self) -> None:
        """
        Initializes various matrices required for calculating the layout of a molecular 
        structure.

        This method computes the initial positions of atoms based on the center of the molecule 
        and bond length. It creates matrices to store bond lengths, link stiffnesses, 
        interaction energies,
        and sums of energy along the X and Y axes. The initialization process involves 
        calculating the distance matrix, determining initial atom positions,
        and preparing matrices for further calculations in the Kamada-Kawai algorithm.

        Steps involved:
        1. Compute the distance matrix to understand the connectivity and distances between 
        atoms.
        2. Determine initial positions for atoms based on a circular layout around the 
        molecule's center, considering unpositioned atoms.
        3. Initialize matrices for bond lengths, spring strengths, and interaction energies 
        between atoms.
        4. Calculate the initial energy matrix based on atom positions and bond lengths.
        """
        self.distance_matrix = self.get_subgraph_distance_matrix(self.atoms)
        length = len(self.atoms)
        radius = Polygon.find_polygon_radius(500, length)
        angle = Polygon.get_central_angle(length)
        a: float = 0.0
        for atom in self.atoms:
            if not atom.positioned:
                self.x_positions[atom] = self.center.x + math.cos(a) * radius
                self.y_positions[atom] = self.center.y + math.sin(a) * radius
            else:
                self.x_positions[atom] = atom.position.x
                self.y_positions[atom] = atom.position.y
            self.positioned[atom] = atom.positioned
            a += angle
        for atom_1 in self.atoms:
            self.length_matrix[atom_1] = {}
            self.spring_strengths[atom_1] = {}
            self.energy_matrix[atom_1] = {}
            self.energy_sums_x[atom_1] = None
            self.energy_sums_y[atom_1] = None
            for atom_2 in self.atoms:
                self.length_matrix[atom_1][atom_2] = self.edge_strength * self.distance_matrix[atom_1][atom_2]
                self.spring_strengths[atom_1][atom_2] = self.edge_strength * self.distance_matrix[atom_1][atom_2] ** -2.0
                self.energy_matrix[atom_1][atom_2] = None
        for atom_1 in self.atoms:
            ux = self.x_positions[atom_1]
            uy = self.y_positions[atom_1]
            d_ex = 0.0
            d_ey = 0.0
            for atom_2 in self.atoms:
                if atom_1 == atom_2:
                    continue
                vx = self.x_positions[atom_2]
                vy = self.y_positions[atom_2]
                denom = 1.0 / math.sqrt((ux - vx) ** 2 + (uy - vy) ** 2)
                self.energy_matrix[atom_1][atom_2] = (self.spring_strengths[atom_1][atom_2] * ((ux - vx) - self.length_matrix[atom_1][atom_2] * (ux - vx) * denom),
                                                      self.spring_strengths[atom_1][atom_2] * ((uy - vy) - self.length_matrix[atom_1][atom_2] * (uy - vy) * denom))
                self.energy_matrix[atom_2][atom_1] = self.energy_matrix[atom_1][atom_2]
                d_ex += self.energy_matrix[atom_1][atom_2][0]
                d_ey += self.energy_matrix[atom_1][atom_2][1]
            self.energy_sums_x[atom_1] = d_ex
            self.energy_sums_y[atom_1] = d_ey
  
  
    def get_kk_layout(self) -> None:
        """
        Initiates the iterative process to find the optimal arrangement of atoms in a molecular 
        structure using the Kamada-Kawai algorithm.

        This method performs iterations until the system's energy falls below a threshold value 
        or the maximum number of iterations is reached. At each iteration,
        the `update()` method is called to move atoms with the highest energy, aiming to 
        minimize the overall system energy through gradual adjustments.

        Description:
        The Kamada-Kawai algorithm is employed to iteratively refine the positions of atoms 
        within a molecular structure, seeking a configuration that minimizes the system's 
        energy. This method orchestrates the iterative process,
        adjusting atom positions based on their energy states until a satisfactory layout is 
        achieved or predefined limits are met. It operates by repeatedly identifying atoms with 
        the highest energy contributions
        and adjusting their positions to reduce strain within the molecular structure, thereby 
        optimizing the layout towards a state of lower potential energy.

        Process:
        - Iterations continue until either the system's energy drops below a specified 
        threshold, indicating an acceptable level of stability, or the maximum iteration count 
        is reached, preventing infinite loops.
        - At each iteration, the atom contributing most significantly to the system's energy is 
        identified, and its position is adjusted to decrease overall energy.
        - Inner iterations within each main iteration further refine the position of the most 
        energetic atom, stopping once the change in energy falls below an inner threshold or the 
        maximum number of inner iterations is reached,
        ensuring fine-tuning of atomic positions for optimal placement.
        - After concluding iterations, final positions are assigned to atoms, marking them as 
        positioned and forcing their placement to prevent further adjustments.

        Outcome:
        - Achieves a stable molecular layout where atomic positions are optimized to minimize 
        energy, reflecting the algorithm's goal of balance and stability.
        - Marks atoms as positioned and forcibly placed, indicating completion and preventing 
        further adjustments, ensuring structural integrity.
        """
        iteration = 0
        max_energy = self.max_energy
        while max_energy > self.threshold and self.max_iteration > iteration:
            iteration += 1
            max_energy_atom, max_energy, d_ex, d_ey = self.highest_energy()
            delta = max_energy
            inner_iteration = 0
            while delta > self.inner_threshold and self.max_inner_iteration > inner_iteration:
                inner_iteration += 1
                self.update(max_energy_atom, d_ex, d_ey)
                delta, d_ex, d_ey = self.energy(max_energy_atom)
        for atom in self.atoms:
            atom.position.x = self.x_positions[atom]
            atom.position.y = self.y_positions[atom]
            atom.positioned = True
            atom.force_positioned = True
   
   
    def energy(self, atom: 'AtomProperties') -> List[float]:
        """
        Calculates the energy of the system for a given atom.

        The energy is defined as the sum of the squares of the energy components along the X and 
        Y axes, as well as the components themselves. This allows for an evaluation of the 
        atom's overall state within the system and its contribution to the total energy.

        Parameters:
        :param atom AtomProperties: 
            The atom for which the energy is calculated.

        Returns List[float]: 
        A list containing:
            - Total energy (sum of the squares of the components),
            - Energy along the X-axis,
            - Energy along the Y-axis.
        """
        energy: List[float] = [self.energy_sums_x[atom]**2 + self.energy_sums_y[atom]**2, \
                               self.energy_sums_x[atom], self.energy_sums_y[atom]]
        return energy
  
  
    def highest_energy(self) -> Tuple['AtomProperties', float, float, float]:
        """
        Identifies the atom with the highest energy among those not yet positioned.

        This method scans through all atoms in the molecular structure to find the one with the 
        greatest energy contribution that has not been positioned yet. It is crucial for the 
        iterative process of optimizing the layout according to the Kamada-Kawai algorithm, as 
        it targets atoms requiring adjustment to minimize overall system energy.

        Returns Tuple[AtomProperties, float, float, float]: 
        A tuple containing:
            - AtomProperties: The atom identified as having the highest energy among those not yet positioned.
            - float: The maximum energy value associated with this atom.
            - float: The energy component along the X-axis for the atom with the highest energy.
            - float: The energy component along the Y-axis for the atom with the highest energy.
        """
        max_energy = 0.0
        max_energy_atom = None
        max_d_ex = 0.0
        max_d_ey = 0.0
        for atom in self.atoms:
            delta, d_ex, d_ey = self.energy(atom)
            if delta > max_energy and not self.positioned[atom]:
                max_energy = delta
                max_energy_atom = atom
                max_d_ex = d_ex
                max_d_ey = d_ey
        return max_energy_atom, max_energy, max_d_ex, max_d_ey
  
  
    def update(self, atom: 'AtomProperties', d_ex: float, d_ey: float) -> None:
        """
        Updates the position of a specified atom based on its energy.

        Parameters:
        :param atom AtomProperties: 
            The atom whose position needs to be updated.
        :param d_ex float: 
            Energy along the X-axis for the atom.
        :param d_ey float: 
            Energy along the Y-axis for the atom.

        Description:
        This method recalculates and adjusts the position of a given atom within the molecular structure based on its current energy state, aiming to minimize the overall system energy through iterative refinement. It incorporates the Kamada-Kawai algorithm principles, adjusting atomic positions to reduce strain and achieve a stable configuration. By considering the energies along the X and Y axes, the method computes new coordinates that reflect a balance between the atom's interactions with other atoms, effectively reducing its contribution to the system's total energy. The process involves calculating forces acting on the atom due to its connections, represented by springs with specific strengths, and updating its position accordingly. The adjustments are made in both X and Y directions, aiming to move the atom towards a state of lower potential energy, thereby contributing to the optimization of the entire molecular layout.
        """
        dxx = 0.0
        dyy = 0.0
        dxy = 0.0
        ux = self.x_positions[atom]
        uy = self.y_positions[atom]
        lengths_array = self.length_matrix[atom]
        strengths_array = self.spring_strengths[atom]
        for atom_2 in self.atoms:
            if atom == atom_2:
                continue
            vx = self.x_positions[atom_2]
            vy = self.y_positions[atom_2]
            length = lengths_array[atom_2]
            strength = strengths_array[atom_2]
            squared_xdiff = (ux - vx) ** 2
            squared_ydiff = (uy - vy) ** 2
            denom = 1.0 / (squared_xdiff + squared_ydiff) ** 1.5
            dxx += strength * (1 - length * squared_ydiff * denom)
            dyy += strength * (1 - length * squared_xdiff * denom)
            dxy += strength * (length * (ux - vx) * (uy - vy) * denom)
        if dxx == 0:
            dxx = 0.1
        if dyy == 0:
            dyy = 0.1
        if dxy == 0:
            dxy = 0.1
        dy = (d_ex / dxx + d_ey / dxy) / (dxy / dxx - dyy / dxy)
        dx = -(dxy * dy + d_ex) / dxx
        self.x_positions[atom] += dx
        self.y_positions[atom] += dy
        d_ex = 0.0
        d_ey = 0.0
        ux = self.x_positions[atom]
        uy = self.y_positions[atom]
        for atom_2 in self.atoms:
            if atom == atom_2:
                continue
            vx = self.x_positions[atom_2]
            vy = self.y_positions[atom_2]
            previous_ex = self.energy_matrix[atom][atom_2][0]
            previous_ey = self.energy_matrix[atom][atom_2][1]
            denom = 1.0 / math.sqrt((ux - vx) ** 2 + (uy - vy) ** 2)
            dx = strengths_array[atom_2] * ((ux - vx) - lengths_array[atom_2] * (ux - vx) * denom)
            dy = strengths_array[atom_2] * ((uy - vy) - lengths_array[atom_2] * (uy - vy) * denom)
            self.energy_matrix[atom][atom_2] = [dx, dy]
            d_ex += dx
            d_ey += dy
            self.energy_sums_x[atom_2] += dx - previous_ex
            self.energy_sums_y[atom_2] += dy - previous_ey
        self.energy_sums_x[atom] = d_ex
        self.energy_sums_y[atom] = d_ey
  
  
    def get_subgraph_distance_matrix(self, atoms: List['AtomProperties']) \
                                            -> Dict[int, Dict[int, float]]:
        """
        Computes the distance matrix between atoms in a subgraph.

        Parameters:
        :param atoms List[AtomProperties]: 
            Specifies the subset of atoms to analyze, focusing the calculation on relevant 
            components of the molecular structure.

        Returns Dict[int, Dict[int, float]]: 
            Provides a comprehensive view of atomic distances, where keys represent atoms and 
            values are dictionaries mapping to other atoms with corresponding shortest path 
            distances. This structured output facilitates targeted adjustments, ensuring that 
            atomic placements minimize overall system energy.

        Description:
        This method calculates the shortest path distances between all pairs of atoms within a 
        given subset of a molecular structure, forming a subgraph. It initializes the distance 
        matrix with infinite distances for all atom pairs except those directly connected, which 
        are set to 1, indicating a bond exists. It then applies a variation of the 
        Floyd-Warshall algorithm to find the shortest paths between all atoms, updating the 
        matrix to reflect the shortest distances found. This process is crucial for 
        understanding the connectivity and spatial relationships within the molecular structure, 
        aiding in layout optimization by identifying the most efficient paths between atoms. The 
        resulting matrix provides insights into the molecular topology, guiding the arrangement 
        of atoms to minimize overall energy and enhance structural stability.
        """
        distance_matrix: Dict[int, Dict[int, float]] = {}
        for atom_1 in atoms:
            if atom_1 not in distance_matrix:
                distance_matrix[atom_1] = {}

            for atom_2 in atoms:
                if self.structure.bond_lookup(atom_1, atom_2):
                    distance_matrix[atom_1][atom_2] = 1
                else:
                    distance_matrix[atom_1][atom_2] = float('inf')

        for atom_1 in atoms:
            for atom_2 in atoms:
                for atom_3 in atoms:
                    if distance_matrix[atom_2][atom_3] > distance_matrix[atom_2][atom_1] + distance_matrix[atom_1][atom_3]:
                        distance_matrix[atom_2][atom_3] = distance_matrix[atom_2][atom_1] + distance_matrix[atom_1][atom_3]
        return distance_matrix