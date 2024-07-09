from collections import defaultdict, OrderedDict
import numpy as np

from ase import Atoms
from ase.data import covalent_radii
from ase.io import write

import matid.geometry
from matid.data import constants

import networkx as nx


class LinkedUnitCollection(dict):
    """Represents a collection of similar cells that are connected in 3D space
    to form a structure, e.g. a surface.

    Essentially this is a special flavor of a regular dictionary: the keys can
    only be a sequence of three integers, and the values should be LinkedUnits.
    """

    def __init__(
        self,
        system,
        cell,
        is_2d,
        dist_matrix_radii_pbc,
        delaunay_threshold=constants.DELAUNAY_THRESHOLD,
        chem_similarity_threshold=constants.CHEM_SIMILARITY_THRESHOLD,
        bond_threshold=constants.BOND_THRESHOLD,
    ):
        """
        Args:
            system(ase.Atoms): A reference to the system from which this
                LinkedUniCollection is gathered.
            cell(ase.Atoms): The prototype cell that is used in finding this
                region.
            is_2d(boolean): Whether this system represents a 2D-material or not.
            delaunay_threshold(float): The maximum allowed size of a tetrahedra
                edge in the Delaunay triangulation of the region..
        """
        self.system = system
        self.cell = cell
        self.is_2d = is_2d
        self.delaunay_threshold = delaunay_threshold
        self.chem_similarity_threshold = chem_similarity_threshold
        self.bond_threshold = bond_threshold
        self.dist_matrix_radii_pbc = dist_matrix_radii_pbc
        self._search_graph = nx.MultiDiGraph()
        self._wrapped_moves = []
        self._index_cell_map = {}
        self._used_points = set()
        self._decomposition = None
        self._inside_indices = None
        self._outside_indices = None
        self._adsorbates = None
        self._substitutions = None
        self._vacancies = None
        self._clusters = None
        self._basis_indices = None
        self._basis_environments = None
        self._translations = None
        self._pos_tol = None
        dict.__init__(self)

    def __setitem__(self, key, value):
        # Transform key to tuple, check length
        try:
            key = tuple(key)
        except Exception:
            raise TypeError(
                "Could not transform the given key '{}' into tuple.".format(key)
            )
        if len(key) != 3:
            raise ValueError(
                "The given coordinate '{}' does not have three components.".format(key)
            )

        # Check that old unit is not overwritten
        if key in dict.keys(self):
            raise ValueError("Overriding existing units is not supported.")

        dict.__setitem__(self, key, value)

    def recreate_valid(self):
        """Used to recreate a new Atoms object, where each atom is created from
        a single unit cell. Atoms that were found not to belong to the periodic
        unit cell are not included.
        """
        recreated_system = Atoms(
            cell=self.system.get_cell(),
            pbc=self.system.get_pbc(),
        )
        for unit in self.values():
            i_valid_indices = np.array([x for x in unit.basis_indices if x is not None])
            if len(i_valid_indices) != 0:
                i_atoms = self.system[i_valid_indices]
                recreated_system += i_atoms

        return recreated_system

    # def get_basis_atom_neighbourhood(self):
    #     """For each atom in the basis calculates the chemical neighbourhood.
    #     The chemical neighbourhood consists of a list of atomic numbers that
    #     are closer than a certain threshold value when the covalent radii is
    #     taken into account.

    #     Args:

    #     Returns:
    #     """
    #     if self._basis_environments is None:
    #         # Multiply the system to get the entire neighbourhood.
    #         cell = self.cell
    #         max_radii = covalent_radii[cell.get_atomic_numbers()].max()
    #         cutoff = max_radii + self.bond_threshold
    #         if self.is_2d:
    #             pbc = [True, True, False]
    #         else:
    #             pbc = [True, True, True]
    #         factors = matid.geometry.get_neighbour_cells(cell.get_cell(), cutoff, pbc)
    #         tvecs = np.dot(factors, cell.get_cell())

    #         # Find the factor corresponding to the original cell
    #         for i_factor, factor in enumerate(factors):
    #             if tuple(factor) == (0, 0, 0):
    #                 tvecs_reduced = np.delete(tvecs, i_factor, axis=0)
    #                 break

    #         pos = cell.get_positions()
    #         disp = matid.geometry.get_displacement_tensor(pos, cell.get_cell())

    #         env_list = []
    #         for i in range(len(cell)):
    #             i_env = self.get_chemical_environment(
    #                 cell, i, disp, tvecs, tvecs_reduced
    #             )
    #             env_list.append(i_env)
    #         self._basis_environments = env_list

    #     return self._basis_environments

    # def get_chem_env_translations(self):
    #     """Used to calculate the translations that are used in calculating the
    #     chemical enviroments.
    #     """
    #     if self._translations is None:
    #         cell = self.system.get_cell()
    #         num = self.system.get_atomic_numbers()
    #         max_radii = covalent_radii[num].max()
    #         cutoff = max_radii + self.bond_threshold
    #         factors = matid.geometry.get_neighbour_cells(cell, cutoff, True)
    #         translations = np.dot(factors, cell)

    #         # Find and remove the factor corresponding to the original cell
    #         for i_factor, factor in enumerate(factors):
    #             if tuple(factor) == (0, 0, 0):
    #                 translations_reduced = np.delete(translations, i_factor, axis=0)
    #                 break
    #         self._translations = (translations, translations_reduced)
    #     return self._translations

    def get_basis_indices(self):
        """Returns the indices of the atoms that were found to belong to a unit
        cell basis in the LinkedUnits in this collection as a single list.

        Returns:
            np.ndarray: Indices of the atoms in the original system that belong
            to this collection of LinkedUnits.
        """
        if self._basis_indices is None:
            # The chemical similarity check is completely skipped if threshold is zero
            # if self.chem_similarity_threshold == 0:
            indices = set()
            for unit in self.values():
                for index in unit.basis_indices:
                    if index is not None:
                        indices.add(index)
            self._basis_indices = indices
            # else:
            # translations, translations_reduced = self.get_chem_env_translations()

            # For each atom in the basis check the chemical environment
            # neighbour_map = self.get_basis_atom_neighbourhood()

            # indices = set()
            # for unit in self.values():
            #     indices |= unit.basis_indices
            # self._basis_indices = np.array(list(indices))

            # Compare the chemical environment near this atom to the one
            # that is present in the prototype cell. If these
            # neighbourhoods are too different, then the atom is not
            # counted as being a part of the region.
            # for i_index, index in enumerate(unit.basis_indices):
            #     if index is not None:
            #         real_environment = self.get_chemical_environment(
            #             self.system,
            #             index,
            #             self.disp_tensor_finite,
            #             translations,
            #             translations_reduced,
            #         )
            #         ideal_environment = neighbour_map[i_index]
            #         chem_similarity = self.get_chemical_similarity(
            #             ideal_environment, real_environment
            #         )
            #         if chem_similarity >= self.chem_similarity_threshold:
            #             indices.add(index)

            # Ensure that all the basis atoms belong to the same cluster.
            # clusters = self.get_clusters()
            # self._basis_indices = np.array(list(indices))

        return self._basis_indices

    # def get_chemical_environment(
    #     self, system, index, disp_tensor_finite, translations, translations_reduced
    # ):
    #     """Get the chemical environment around an atom. The chemical
    #     environment is quantified simply by the number of different species
    #     around a certain distance when the covalent radii have been considered.
    #     """
    #     # Multiply the system to get the entire neighbourhood.
    #     num = system.get_atomic_numbers()
    #     n_atoms = len(system)
    #     seed_num = num[index]

    #     neighbours = defaultdict(lambda: 0)
    #     for j in range(n_atoms):
    #         j_num = num[j]
    #         ij_disp = disp_tensor_finite[index, j, :]

    #         if index == j:
    #             trans = translations_reduced
    #         else:
    #             trans = translations

    #         D_trans = trans + ij_disp
    #         D_trans_len = np.linalg.norm(D_trans, axis=1)
    #         ij_radii = covalent_radii[seed_num] + covalent_radii[j_num]
    #         ij_n_neigh = np.sum(D_trans_len - ij_radii <= self.bond_threshold)

    #         neighbours[j_num] += ij_n_neigh

    #     return neighbours

    # def get_chemical_similarity(self, ideal_env, real_env):
    #     """Returns a metric that quantifies the similarity between two chemical
    #     environments. Here the metric is defined simply by the number
    #     neighbours that are found to be same as in the ideal environmen within
    #     a certain radius.
    #     """
    #     max_score = sum(ideal_env.values())

    #     score = 0
    #     for ideal_key, ideal_value in ideal_env.items():
    #         real_value = real_env.get(ideal_key)
    #         if real_value is not None:
    #             score += min(real_value, ideal_value)

    #     return score / max_score

    def get_all_indices(self):
        """Get all the indices that are present in the full system."""
        return set(range(len(self.system)))

    def get_connected_directions(self):
        """During the tracking of the region the information about searches
        that matched an atom twice but with a negated multiplier are stored.
        """
        G = self._search_graph
        dir_vectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        directions = set([0, 1, 2])
        for node in G.nodes():
            node_edges = G.in_edges(node, data=True)
            dir_to_remove = set()
            for direction in directions:
                dir_vector = dir_vectors[direction]
                positive = False
                negative = False
                for edge in node_edges:
                    multiplier = edge[2]["multiplier"]
                    if np.array_equal(multiplier, dir_vector):
                        positive = True
                    if np.array_equal(multiplier, -dir_vector):
                        negative = True
                    if positive and negative:
                        break
                if positive and negative:
                    dir_to_remove.add(direction)
            directions -= dir_to_remove

        connected_directions = np.array([True, True, True])
        connected_directions[list(directions)] = False
        return connected_directions


class LinkedUnit:
    """Represents a cell that is connected to others in 3D space to form a
    structure, e.g. a surface.
    """

    def __init__(
        self,
        index,
        seed_index,
        seed_coordinate,
        cell,
        basis_indices,
        substitutions,
        vacancies,
    ):
        """
        Args:
            index(tuple of three ints):
            seed_index(int):
            seed_coordinate():
            cell(np.ndarray): Cell for this unit. Can change from unit to unit.
            basis_indices(sequence of ints and Nones): A sequence where there
                is an index or None for each atom that is supposed to be in the
                basis
            substitute_indices(sequence of ints and Nones): If basis atom is
                replaced by a foreign atom, the index of the substitutional atom is
                here.
        """
        self.index = index
        self.seed_index = seed_index
        self.seed_coordinate = seed_coordinate
        self.cell = cell
        self.basis_indices = basis_indices
        self.substitutions = substitutions
        self.vacancies = vacancies


class Substitution:
    """Represents a substitutional point defect."""

    def __init__(self, index, position, original_element, substitutional_element):
        """
        Args:
            index (int): Index of the subtitutional defect in the original
                system.
            original (ase.Atom): The original atom that was substituted
            substitution (ase.Atom): The substitutional atom

        """
        self.index = index
        self.original_element = original_element
        self.substitutional_element = substitutional_element
