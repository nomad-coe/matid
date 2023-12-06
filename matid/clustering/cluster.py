import numpy as np
from ase import Atoms

import matid.geometry


class Cluster:
    """
    Represents a part of a bigger system.
    """

    def __init__(
        self,
        indices=None,
        species=None,
        region=None,
        dimensionality=None,
        cell=None,
        system=None,
        distances=None,
        radii=None,
        bond_threshold=None,
    ):
        if isinstance(indices, list):
            self.indices = indices
        else:
            self.indices = list(indices)
        self.species = species

        self._region = region
        self._dimensionality = dimensionality
        self._cell = cell
        self._system = system
        self._distances = distances
        self._radii = radii
        self._merged = False
        self._bond_threshold = bond_threshold
        self._distance_matrix_radii_mic = None

    def __len__(self):
        return len(self.indices)

    def _get_distance_matrix_radii_mic(self) -> np.ndarray:
        """Retrieves the distance matrix with subtracted radii for this cluster."""
        if self._distance_matrix_radii_mic is None:
            self._distance_matrix_radii_mic = self._distances.dist_matrix_radii_mic[np.ix_(self.indices, self.indices)]
        return self._distance_matrix_radii_mic

    def get_cell(self) -> Atoms:
        """Used to fetch the prototypical cell for this cluster if one exists."""
        if self._cell:
            return self._cell
        if self._region:
            return self._region.cell
        return None

    def get_atoms(self) -> Atoms:
        """Returns the ase.Atoms object for this cluster."""
        return self._system[self.indices]

    def get_dimensionality(self) -> int:
        """Shortcut for fetching the dimensionality of the cluster using
        matid.geometry.get_dimensionality and the radii + bond thresholds that
        were used during the clustering.
        """
        if self._dimensionality is None:
            self._dimensionality = matid.geometry.get_dimensionality(
                self.get_atoms(),
                self._bond_threshold,
                dist_matrix_radii_mic_1x=self._get_distance_matrix_radii_mic(),
            )
        return self._dimensionality