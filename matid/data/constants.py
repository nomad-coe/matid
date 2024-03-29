# The variable SYMMETRY_TOL controls the precision used by spglib in order to
# find symmetries
SYMMETRY_TOL = 0.4  # unit: angstrom

# The threshold for a system to be considered "flat". Used e.g. when
# determining if a 2D structure is purely 2-dimensional to allow extra rigid
# transformations that are improper in 3D but proper in 2D.
FLAT_DIM_THRESHOLD = 0.1

# An ordered list of Wyckoff letters
WYCKOFF_LETTERS = list("abcdefghijklmnopqrstuvwxyzA")
WYCKOFF_LETTER_POSITIONS = {
    letter: positions for positions, letter in enumerate(WYCKOFF_LETTERS)
}

# ===============================================================================
# Constants for classification
MAX_SINGLE_CELL_SIZE = 5
MAX_2D_CELL_HEIGHT = 5
MAX_CELL_SIZE = [12]
REL_POS_TOL = [
    0.25,
    0.75,
]  # Position tolerances. Given relative to minimum distance between two atoms in the system.
ANGLE_TOL = 20
CLUSTER_THRESHOLD = 3.5
CRYSTALLINITY_THRESHOLD = 0.25
DELAUNAY_THRESHOLD = 1.5
BOND_THRESHOLD = 0.75
CHEM_SIMILARITY_THRESHOLD = 0.4
CELL_SIZE_TOL = 0.25
MIN_COVERAGE = 0.5
