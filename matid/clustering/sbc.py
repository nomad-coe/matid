from collections import defaultdict

import numpy as np
import ase.geometry

import matid.geometry
from matid.clustering.cluster import Cluster
from matid.core.periodicfinder import PeriodicFinder


class SBC:
    """
    Class for performing Symmetry-based clustering (SBC).

    You can apply this class for partitioning a larger material into grains, a
    heterostructure into it's component etc. The clustering is based on finding
    periodically repeating motifs, and as such it is not suitable for e.g.
    finding molecules. Any atoms that do not have enough periodic repetitions
    will be returned as isolated clusters.
    """

    def __init__(self, seed=7):
        """
        Args:
            seed(int): The seed that is used for random number generation.
        """
        self.rng = np.random.default_rng(seed)

    def _merge_clusters(
        self, system, clusters, merge_threshold, distances, bond_threshold
    ):
        """
        Used to merge clusters that have the same species and share atoms.
        """

        def merge(system, a, b):
            """
            Merges the given two clusters.
            """
            # If there are conflicting species in the regions that are merged,
            # only the species from the larger are kept during the merge. This
            # helps getting rid of artifical regions at the interface between
            # two clusters.
            atomic_numbers = system.get_atomic_numbers()
            if len(a.indices) > len(b.indices):
                target = a
                source = b
            else:
                target = b
                source = a
            common = set(
                filter(lambda x: atomic_numbers[x] in target.species, source.indices)
            )
            final_indices = set(target.indices).union(common)

            # TODO: The merged cluster simply inherits the largest of the two
            # regions. There currently is no mechanism for mergin two regions
            # together.
            sorted_regions = sorted(
                [a._region, b._region],
                key=lambda x: -1 if x is None else len(x.get_basis_indices()),
            )
            largest_region = sorted_regions[-1]

            return Cluster(
                final_indices,
                target.species,
                largest_region,
                system=system,
                distances=distances,
                bond_threshold=bond_threshold,
            )

        isolated_clusters = []
        while True:
            if len(clusters) == 0 or clusters[0]._merged:
                break
            i_cluster = clusters.pop(0)
            i_indices = set(i_cluster.indices)

            # Check overlap with all other non-isolated clusters
            isolated = True
            if len(clusters):
                overlaps = [
                    (j, len(i_indices.intersection(set(j_cluster.indices))))
                    for j, j_cluster in enumerate(clusters)
                ]
                overlaps = sorted(overlaps, key=lambda x: x[1], reverse=True)
                best_overlap = overlaps[0][1]
                best_grain = overlaps[0][0]
                target_cluster = clusters[best_grain]

                # Find the biggest overlap and if it is large enough, merge
                # these clusters together. Large overlap indicates that they are
                # actually part of the same component that the larger region
                # represents a better description.
                best_overlap_score = max(
                    best_overlap / len(i_indices),
                    best_overlap / len(target_cluster.indices),
                )

                if best_overlap_score > merge_threshold:
                    merged = merge(system, i_cluster, target_cluster)
                    merged._merged = True
                    isolated = False
                    clusters.pop(best_grain)
                    clusters.append(merged)
            # Component without enough overlap is saved as it is.
            if isolated:
                isolated_clusters.append(i_cluster)

        return isolated_clusters + clusters

    def _localize_clusters(self, system, clusters, merge_radius, distances):
        """
        Used to resolve overlaps between clusters by assigning the overlapping
        atoms to the "nearest" cluster.

        Args:
            system (ase.Atoms): The original system.
            clusters (list of Clusters): The clusters to localize.
            merge_radius (float): The radius to consider when assigning atoms
              to a particular cluster.

        Returns:
            List of Clusters that no longer have overlapping atoms as each atom
            has been assigned to the nearest cluster.
        """
        # Get all overlapping atoms, and the regions with which they overlap
        overlap_map = defaultdict(list)
        for i in range(len(system)):
            for cluster in clusters:
                if i in cluster.indices:
                    overlap_map[i].append(cluster)

        # Assign each overlapping atom to the cluster that is "nearest". Notice
        # that we do not update the regions during the process.
        # positions = system.get_positions()
        for i, i_clusters in overlap_map.items():
            if len(i_clusters) > 1:
                surrounding_indices = set(
                    np.argwhere(distances.dist_matrix_radii_mic[i, :] < merge_radius)[
                        :, 0
                    ]
                )
                max_near = 0
                max_cluster = i_clusters[0]
                for cluster in i_clusters:
                    n_near = len(set(cluster.indices).intersection(surrounding_indices))
                    if n_near > max_near:
                        max_near = n_near
                        max_cluster = cluster

                for cluster in i_clusters:
                    ind_set = set(cluster.indices)
                    if cluster != max_cluster:
                        ind_set.remove(i)
                    cluster.indices = list(ind_set)
        return clusters

    def _clean_clusters(self, clusters, bond_threshold):
        """
        Cleans out any "dangling" atoms from the cluster by examining chemical
        connectivity. Required because the periodic search does not care about
        the chemical connectivity of all cells in the region, only the prototype
        cell is checked.

        Args:
            system (ase.Atoms): The original system.
            clusters (list of Clusters): The clusters to localize.

        Returns:
            List of Clusters where any outlier atoms have been removed.
        """
        for cluster in clusters:
            dbscan_clusters = matid.geometry.get_clusters(
                cluster._get_distance_matrix_radii_mic(),
                bond_threshold,
                min_samples=1,
            )
            largest_indices = max(dbscan_clusters, key=lambda x: len(x))
            cluster.indices = np.array(cluster.indices)[largest_indices].tolist()
        return clusters

    def get_clusters(
        self,
        system,
        angle_tol=20,
        max_cell_size=6,
        pos_tol=0.7,
        merge_threshold=0.5,
        merge_radius=1,
        bond_threshold=0.65,
        overlap_threshold=-0.1,
        radii="covalent",
    ):
        """
        Used to detect and return structurally separate clusters within the
        given system.

        Args:
            system (ase.Atoms): The structure to partition.
            angle_tol (float): angle_tol parameter for PeriodicFinder
            max_cell_size (float): max_cell_size parameter for PeriodicFinder.get_region
            pos_tol (float): pos_tol parameter for PeriodicFinder.get_region
            merge_threshold (float): A threshold for merging two clusters
                together. Give as a fraction of shared atoms. Value of 1 would
                mean that clusters are never merged, value of 0 means that they
                are merged always when at least one atom is shared.
            merge_radius (float): Radius for finding nearby atoms when deciding
                which cluster is closest. The atomic radii are subtracted from
                distances. Given in angstroms.
            bond_threshold(float): Used to control the connectivity threshold
                for defining a chemical connection between atoms. Controls e.g.
                what types of unit cells are accepted and how outliers are
                removed from clusters.
            overlap_threshold(float): Used to exclude non-physical cells by
                checking overlap of atoms. Overlap between two atoms is
                calculated by subtracting atomic radii from the distance between
                the atoms.
            radii(str|np.ndarray): The radii to use for atoms. Use either a preset
                or a custom list of atomic radii for each atom. The available presets are:

                    - covalent: Covalent radii from DOI:10.1039/B801115J
                    - vdw: van Der Waals radii from DOI:10.1039/C3DT50599E
                    - vdw_covalent: preferably van Der Waals radii, covalent if vdw
                    not defined.

        Returns:
            A list of Clusters.
        """
        # Copy the system to avoid mutating the original
        system_copy = system.copy()

        # For each lattice basis vector that has a zero length, and which is not
        # periodic, we extend the basis vectors to fully contain all the atoms.
        # This is to ensure that the code that requires scaled positions works
        # correctly.
        pbc = system_copy.get_pbc()
        basis = system_copy.get_cell()
        requires_completion = False
        for i in range(3):
            if not basis[i, :].any():
                if not pbc[i]:
                    requires_completion = True
                else:
                    raise ValueError(
                        "Cannot process system with zero-volume cell and periodic boundaries."
                    )
        if requires_completion:
            system_copy.set_cell(ase.geometry.complete_cell(basis))
            scaled_positions = system_copy.get_scaled_positions()
            new_cell = system_copy.get_cell()
            for i in range(3):
                if not pbc[i]:
                    max_pos = scaled_positions[:, i].max()
                    min_pos = scaled_positions[:, i].min()
                    new_cell[i, :] *= (max_pos - min_pos) + 1
            system_copy.set_cell(new_cell)
            system_copy.center()

        # Positions are wrapped
        system_copy.wrap()

        atomic_numbers = system.get_atomic_numbers()
        radii = matid.geometry.get_radii(radii, atomic_numbers)

        # Calculate the distances here once if they have not been provided.
        distances = matid.geometry.get_distances(system_copy, radii)

        # Iteratively search for new clusters until whole system is covered
        periodic_finder = PeriodicFinder(
            angle_tol=angle_tol, chem_similarity_threshold=0
        )
        indices = set(list(range(len(system_copy))))
        clusters = []
        while len(indices) != 0:
            i_seed = self.rng.choice(list(indices), 1)[0]
            i_grain, mask = periodic_finder.get_region(
                system_copy,
                seed_index=i_seed,
                max_cell_size=max_cell_size,
                pos_tol=pos_tol,
                bond_threshold=bond_threshold,
                overlap_threshold=overlap_threshold,
                distances=distances,
                return_mask=True,
            )

            # All neighbours that the periodic finder has tested are removed
            # from the search. This significantly helps with the scaling of the
            # clustering.
            tested_indices = set(np.arange(len(mask))[mask])
            indices -= tested_indices

            # If a grain is found, it is added as a single cluster and removed
            # from the search
            if i_grain is not None:
                i_indices = {i_seed}
                i_indices.update(i_grain.get_basis_indices())
                i_species = set(atomic_numbers[list(i_indices)])
                clusters.append(
                    Cluster(
                        i_indices,
                        i_species,
                        i_grain,
                        system=system_copy,
                        distances=distances,
                        radii=radii,
                        bond_threshold=bond_threshold,
                    )
                )
                indices -= i_indices

        # Check overlaps of the regions. For large overlaps the grains are
        # merged (the real region was probably cut into pieces by unfortunate
        # selection of the seed atom)
        clusters = self._merge_clusters(
            system_copy, clusters, merge_threshold, distances, bond_threshold
        )

        # Any remaining overlaps are resolved by assigning atoms to the
        # "nearest" cluster
        clusters = self._localize_clusters(
            system_copy, clusters, merge_radius, distances
        )

        # Any atoms that are not chemically connected to the region will be
        # excluded.
        clusters = self._clean_clusters(clusters, bond_threshold)

        return clusters
