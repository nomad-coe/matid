import numpy as np
from systax.data.element_data import get_symbols, get_numbers
from systax.core.lattice import Lattice


class System(object):
    """Represents an atomic system.

    Args:
        lattice (3x3 ndarray): The lattice for this system.
        positions (ndarray): The relative or cartesian positions of the atoms.
            The type is controlled by 'coords_are_cartesian'.
        species (ndarray): The atomic numbers or symbols of the atoms.
        charges (ndarray): The charges of the atoms.
        coords_are_cartesian (bool): Whether the given coordinates are
            cartesian. If false, the coordinates are taken to be relative to the
            cell vectors.
        wyckoff_letters (ndarray): The Wyckoff letters for the atoms.
        equivalent atoms (ndarray): A list that contains as integer for each
            atom in the system. Same integer for different atoms means that they
            are symmetrically equivalent.
    """
    def __init__(
            self,
            positions,
            species,
            lattice=None,
            charges=None,
            coords_are_cartesian=False,
            wyckoff_letters=None,
            equivalent_atoms=None,
            periodicity=None
            ):

        self._numbers = None
        self._symbols = None
        self._cartesian_pos = None
        self._relative_pos = None
        self.periodicity = periodicity
        positions = np.array(positions)
        if coords_are_cartesian:
            self._cartesian_pos = positions
        else:
            self._relative_pos = positions

        species = np.array(species)
        if isinstance(species[0], np.str):
            self._symbols = species
        elif isinstance(species.item(0), (int, float)):
            self._numbers = species
        else:
            raise ValueError("The type is {}".format(type(species[0])))

        if isinstance(lattice, Lattice):
            self.lattice = Lattice(lattice.matrix)
        if lattice is None:
            self.lattice = None
        else:
            self.lattice = Lattice(lattice)

        if charges is not None:
            charges = np.array(charges)

        self.system_info = None
        self.charges = charges
        self.wyckoff_letters = wyckoff_letters
        self.equivalent_atoms = equivalent_atoms

        self._displacement_tensor = None
        self._distance_matrix = None
        self._inverse_distance_matrix = None

    @classmethod
    def fromatoms(cls, atoms):
        """Create a System object from an ASE Atoms object. This allows users
        to use a familiar library with support for loading configurations from
        all sorts of different files (ase.io), but still enables some
        performance tweaks that are present in the System object, like caching.

        Args:
            atoms (ase.Atoms): The Atoms object from which to create the system.

        Returns:
            System: The new System object identical to the given Atoms.
        """
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        atomic_numbers = atoms.get_atomic_numbers()
        pbc = atoms.get_pbc()

        system = System(
            positions=positions,
            species=atomic_numbers,
            lattice=cell,
            charges=atoms.get_initial_charges(),
            coords_are_cartesian=True,
            periodicity=pbc,
        )

        return system

    @property
    def relative_pos(self):
        if self._relative_pos is None:
            self._relative_pos = np.linalg.solve(self.lattice._matrix.T, self._cartesian_pos.T).T
        return self._relative_pos

    @property
    def cartesian_pos(self):
        if self._cartesian_pos is None:
            self._cartesian_pos = self._relative_pos.dot(self.lattice._matrix.T)
        return self._cartesian_pos

    @property
    def symbols(self):
        if self._symbols is None:
            self._symbols = get_symbols(self._numbers)
        return self._symbols

    @property
    def charge(self):
        return np.sum(self.charges)

    @property
    def numbers(self):
        if self._numbers is None:
            self._numbers = get_numbers(self._symbols)
        return self._numbers

    @property
    def displacement_tensor(self):
        """A matrix where the entry A[i, j, :] is the vector
        self.cartesian_pos[i] - self.cartesian_pos[j]
        """
        if self._displacement_tensor is None:
            pos = self.cartesian_pos

            # Add new axes so that broadcasting works nicely
            disp_tensor = pos[:, None, :] - pos[None, :, :]

            self._displacement_tensor = disp_tensor
        return self._displacement_tensor

    @property
    def distance_matrix(self):
        if self._distance_matrix is None:
            displacement_tensor = self.displacement_tensor
            distance_matrix = np.linalg.norm(displacement_tensor, axis=2)
            self._distance_matrix = distance_matrix
        return self._distance_matrix

    @property
    def inverse_distance_matrix(self):
        if self._inverse_distance_matrix is None:
            distance_matrix = self.distance_matrix
            with np.errstate(divide='ignore'):
                inv_distance_matrix = np.reciprocal(distance_matrix)
            self._inverse_distance_matrix = inv_distance_matrix
        return self._inverse_distance_matrix

    def wrap_positions(self, precision=1E-5):
        """Wrap the relative positions so that each element in the array is within the
        half-closed interval [0, 1)

        By wrapping values near 1 to 0 we will have a consistent way of
        presenting systems.
        """
        self._relative_pos %= 1

        abs_zero = np.absolute(self._relative_pos)
        abs_unity = np.absolute(abs_zero-1)

        near_zero = np.where(abs_zero < precision)
        near_unity = np.where(abs_unity < precision)

        self._relative_pos[near_unity] = 0
        self._relative_pos[near_zero] = 0

        self._reset_cartesian_pos

    def translate(self, translation, relative=True):
        """Translates the relative positions by the given translation.

        Args:
            translation (1x3 numpy.array): The translation to apply.
            relative (bool): True if given translation is relative to cell
                vectors.
        """
        if relative:
            self._relative_pos += translation
            self._reset_cartesian_pos()
        else:
            self._cartesian_pos += translation
            self._reset_relative_pos()

    def _reset_relative_pos(self):
        self._relative_pos = None
        self._distance_matrix = None
        self._displacement_tensor = None
        self._inverse_distance_matrix = None

    def _reset_cartesian_pos(self):
        self._cartesian_pos = None
        self._distance_matrix = None
        self._displacement_tensor = None
        self._inverse_distance_matrix = None

    def __len__(self):
        return len(self.numbers)
