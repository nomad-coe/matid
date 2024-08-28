---
sidebar_position: 4
sidebar_label: Geometry Module
---
# Geometry Module
This module defines functions for deriving geometry related quantities from
a atomic system.

## get\_dimensionality

```python
def get_dimensionality(system,
                       cluster_threshold=CLUSTER_THRESHOLD,
                       dist_matrix_radii_mic_1x=None,
                       return_clusters=False,
                       radii="covalent")
```

Used to calculate the dimensionality of a system with a modified
Topology Scaling Algorithm (TSA) (Michael Ashton, Joshua Paul, Susan B.
Sinnott, and Richard G. Hennig Phys. Rev. Lett. 118, 106101).

**Arguments**:

- `system` _ASE.Atoms_ - The system for which the dimensionality is
  evaluated.
- `cluster_threshold(float)` - The epsilon value for the DBSCAN algorithm
  that is used to identify clusters within the unit cell.
- `dist_matrix_radii_pbc` _np.ndarray_ - A precalculated distance matrix
  that takes in to account the periodicity and has the covalent radii of
  the atoms already subtracted.
- `return_clusters(boolean)` - Whether the clusters are returned
- `radii(str|np.ndarray)` - The radii to use for atoms. Use either a preset
  or a custom list or atomic radii where the atomic number is used as an
  index. The available presets are:
  
  - covalent: Covalent radii from DOI:10.1039/B801115J
  - vdw: van Der Waals radii from DOI:10.1039/C3DT50599E
  - vdw_covalent: preferably van Der Waals radii, covalent if vdw
  not defined.
  

**Returns**:

- `int|none` - The dimensionality of the system. If the dimensionality can't be
  evaluated because the system has multiple disconnected components in
  the original cell, None is returned.
- `list` - A list of clusters. Each entry in the list contains the indices
  of atoms in a cluster.

## get\_moments\_of\_inertia

```python
def get_moments_of_inertia(system, weight=True)
```

Calculates geometric inertia tensor, i.e., inertia tensor but with
all masses are set to 1.

I_ij = sum_k m_k (delta_ij * r_k^2 - x_ki * x_kj)
with r_k^2 = x_k1^2 + x_k2^2 x_k3^2

**Arguments**:

  system(ASE Atoms): Atomic system.
  

**Returns**:

  (np.ndarray, np.ndarray): The eigenvalues and eigenvectors of the
  geometric inertia tensor.

## get\_center\_of\_mass

```python
def get_center_of_mass(system)
```

Calculates the center of mass and also takes the periodicity of the
system into account.

The algorithm is replicated from:
https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions

**Arguments**:

- `system(ase.Atoms)` - The system for which the center of mass is
  calculated.
  

**Returns**:

- `np.ndarray` - The cartesian positions of the center of mass in the given
  system.

## get\_space\_filling

```python
def get_space_filling(system)
```

Calculates the ratio of vacuum to filled space by assuming covalent
radii for the atoms.

**Arguments**:

- `system(ASE.Atoms)` - Atomic system.
  

**Returns**:

- `float` - The ratio of occupied volume to the cell volume.

## make\_random\_displacement

```python
def make_random_displacement(system, delta, rng=None)
```

Dislocate every atom in the given system in a random direction but by
the same amount. Uses an internal random number generator to avoid touching
the global numpy.random.seed()-function.

**Arguments**:

- `system(ASE.Atoms)` - The system for which the displacement are performed.
- `delta(float)` - The magnitude of the displacements.
- `rng(np.random.RandomState)` - Random number generator.

## get\_extended\_system

```python
def get_extended_system(system, cutoff=0)
```

Given a system and a cutoff value, returns a new system which has been
extended so that for each atom the neighbourhood within the cutoff radius is
present, taking periodic boundary conditions into account.

**Arguments**:

- `system(ASE.Atoms)` - System to extend
- `cutoff(float)` - Radial cutoff
  

**Returns**:

  ExtendedSystem object.

## get\_clusters

```python
def get_clusters(dist_matrix, threshold, min_samples=1)
```

Used to detect clusters with the DBSCAN algorithm.

**Arguments**:

- `dist_matrix(np.ndarray)` - The 2D distance matrix from which the clusters
  are calculated.
- `threshold(float)` - The epsilon threshold value for the clustering
- `min_samples(int)` - The minimum allowed cluster size.
  

**Returns**:

- `list` - A list of clusters, where each cluster is a list of indices for
  the elements belonging to the cluster.

## get\_covalent\_distances

```python
def get_covalent_distances(system, mic=True)
```

Returns a distance matrix where the covalent radii have been taken into
account. Clips negative values to be zero.

## get\_biggest\_gap\_indices

```python
def get_biggest_gap_indices(coordinates)
```

Given the list of coordinates for one axis, this function will find the
maximum gap between them and return the index of the bottom and top
coordinates. The bottom and top are defined as:

===       ===============    --->
^top    ^bot               ^axis direction

## get\_dimensions

```python
def get_dimensions(system, vacuum_gaps)
```

Given a system with vacuum gaps, calculate its dimensions in the
directions with vacuum gaps by also taking into account the atomic radii.

## get\_wrapped\_positions

```python
def get_wrapped_positions(scaled_pos, precision=1e-5)
```

Wrap the given relative positions so that each element in the array
is within the half-closed interval [0, 1)

By wrapping values near 1 to 0 we will have a consistent way of
presenting systems.

## find\_mic

```python
def find_mic(D, cell, pbc, max_distance=None)
```

Finds the minimum-image representation of vector(s) D.

**Arguments**:

- `D(np.ndarray)` - A nx3 array of vectors.
- `cell(np.ndarray)` - A valid ase Atoms cell definition.
- `pbc(np.ndarray)` - A 3x1 boolean array for the periodic boundary
  conditions.
- `mic_copies(np.ndarray)` - The maximum number of periodic copies to
  consider in each direction. If not specified, the maximum possible
  number of copies is determined and used.
  

**Returns**:

- `np.ndarray` - The minimum image versions of the given vectors.
- `np.ndarray` - The shifts corresponding to the indices of the neighbouring
  cells in which the vectors were found.

## get\_neighbour\_cells

```python
def get_neighbour_cells(cell, cutoff, pbc)
```

Given a cell and a cutoff, returns the indices of the copies of the cell
which have to be searched in order to reach atom within the cutoff
distance.

The number of neighboring images to search in each direction is equal to
the ceiling of the cutoff distance (defined above) divided by the length of
the projection of the lattice vector onto its corresponding surface normal.
a's surface normal vector is e.g.  b x c / (|b| |c|), so this projection is
(a . (b x c)) / (|b| |c|).  The numerator is just the lattice volume, so
this can be simplified to V / (|b| |c|). This is rewritten as V |a| / (|a|
|b| |c|) for vectorization purposes.

**Arguments**:

  

## get\_mic\_vector

```python
def get_mic_vector(w, v, cell)
```

Used to calculate the minimum image version of a vector in the given
cell.

**Arguments**:

- `rel_vector(np.ndarray)` - Relative vector in the cell basis:
  

**Returns**:

- `np.ndarray` - The MIC version of the given vector
- `np.ndarray` - The shift that corresponds to the neighbouring cell in
  which the copy was found.

## expand\_pbc

```python
def expand_pbc(pbc)
```

Used to expand a pbc definition into an array of three booleans.

**Arguments**:

  pbc(boolean or a list of booleans): The periodicity of the cell. This
  can be any of the values that is also supprted by ASE, namely: a
  boolean or a list of three booleans.
  

**Returns**:

  np.ndarray of booleans: The periodicity expanded as an explicit list of
  three boolean values.

## change\_basis

```python
def change_basis(positions, basis, offset=None)
```

Transform the given cartesian coordinates to a basis that is defined by
the given basis and origin offset.

**Arguments**:

- `positions(np.ndarray)` - Positions in cartesian coordinates.
- `basis(np.ndarray)` - Basis to which to transform.
- `offset(np.ndarray)` - Offset of the origins. A vector from the old basis
  origin to the new basis origin.

**Returns**:

- `np.ndarray` - Relative positions in the new basis

## get\_positions\_within\_basis

```python
def get_positions_within_basis(system,
                               basis,
                               origin,
                               tolerance,
                               mask=[True, True, True],
                               pbc=True)
```

Used to return the indices of positions that are inside a certain basis.
Also takes periodic boundaries into account.

**Arguments**:

- `system(ASE.Atoms)` - System from which the positions are searched.
- `basis(np.ndarray)` - New basis vectors.
- `origin(np.ndarray)` - New origin of the basis in cartesian coordinates.
- `tolerance(float)` - The matching tolerance in angstrom.
  mask(sequence of bool): Mask for selecting the basis's to consider.
  pbc(sequence of bool): The periodicity of the system.
  

**Returns**:

  sequence of int: Indices of the atoms within this cell in the given
  system.
- `np.ndarray` - Relative positions of the found atoms.
- `np.ndarray` - The index of the periodic copy in which the position was
  found.

## get\_matches

```python
def get_matches(system, cell_list, positions, numbers, tolerance)
```

Given a system and a list of cartesian positions and atomic numbers,
returns a list of indices for the atoms corresponding to the given
positions with some tolerance.

**Arguments**:

- `system(ASE.Atoms)` - System where to search the positions
- `cell_list(CellList)` - The cell list for an appropriately extended version
  of the system.
- `positions(np.ndarray)` - Positions to match in the system.
- `tolerance(float)` - Maximum allowed distance for matching.
  

**Returns**:

- `np.ndarray` - indices of matched atoms
- `list` - list of substitutions
- `list` - list of vacancies
- `np.ndarray` - for each searched position, an integer array representing
  the number of the periodic copy where the match was found.

## get\_matches\_simple

```python
def get_matches_simple(system, cell_list, positions, numbers, tolerance)
```

Given a system and a list of cartesian positions and atomic numbers,
returns a list of indices for the atoms corresponding to the given
positions with some tolerance.

**Arguments**:

- `system(ASE.Atoms)` - System where to search the positions
- `cell_list(CellList)` - The cell list for an appropriately extended version
  of the system.
- `positions(np.ndarray)` - Positions to match in the system.
- `tolerance(float)` - Maximum allowed distance for matching.

**Returns**:

- `list` - list of matched atoms or None is nothing was matched.

## get\_cell\_list

```python
def get_cell_list(positions, cell, pbc, extension, cutoff)
```

Given a system and a cutoff value, returns a cell list object.

**Arguments**:

- `positions(np.ndarray)` - Cartesian positions
- `cell(np.ndarray)` - Cell as 3x3 array
- `pbc(np.ndarray)` - Periodic boundary conditions as array of three booleans
- `extension(float)` - How much the system should be extended for the search.
- `cutoff(float)` - Radial cutoff
  

**Returns**:

  CellList object.

## to\_scaled

```python
def to_scaled(cell, positions, wrap=False, pbc=False)
```

Used to transform a set of positions to the basis defined by the
cell of this system.

**Arguments**:

- `system(ASE.Atoms)` - Reference system.
- `positions` _numpy.ndarray_ - The positions to scale
- `wrap` _numpy.ndarray_ - Whether the positions should be wrapped
  inside the cell.
  

**Returns**:

- `numpy.ndarray` - The scaled positions

## to\_cartesian

```python
def to_cartesian(cell, scaled_positions, wrap=False, pbc=False)
```

Used to transform a set of relative positions to the cartesian basis
defined by the cell of this system.

**Arguments**:

- `cell` _numpy.ndarray_ - 3x3 array with lattice vectors as rows.
- `positions` _numpy.ndarray_ - The positions to scale. These postiions
  should have a shape of [n, 3], where n is the number of positions.
- `wrap` _numpy.ndarray_ - Whether the positions should be wrapped
  inside the cell.
  

**Returns**:

- `numpy.ndarray` - The cartesian positions

## translate

```python
def translate(system, translation, relative=False)
```

Translates the positions by the given translation.

**Arguments**:

- `translation` _1x3 numpy.array_ - The translation to apply.
- `relative` _bool_ - True if given translation is relative to cell
  vectors.

## get\_closest\_direction

```python
def get_closest_direction(vec, directions, normalized=False)
```

Used to return the direction that is most parallel to a given one.

**Arguments**:

  

# Intervals Class

Handles list of intervals.

This class allows sorting and adding up of intervals and taking into
account if they overlap.

## \_\_init\_\_

```python
def __init__(intervals=None)
```

**Arguments**:

- `intervals` - List of intervals that are added.

## add\_interval

```python
def add_interval(a, b)
```

Add one interval.

**Arguments**:

  a, b: Start and end of interval. The order does not matter.

## add\_intervals

```python
def add_intervals(intervals)
```

Add list of intervals.

**Arguments**:

- `intervals` - List of intervals that are added.

## set\_intervals

```python
def set_intervals(intervals)
```

Set list of intervals.

**Arguments**:

- `intervals` - List of intervals that are set.

## remove\_interval

```python
def remove_interval(i)
```

Remove one interval.

**Arguments**:

- `i` - Index of interval that is removed.

## get\_intervals

```python
def get_intervals()
```

Returns the intervals.

## get\_intervals\_sorted\_by\_start

```python
def get_intervals_sorted_by_start()
```

Returns list with intervals ordered by their start.

## get\_intervals\_sorted\_by\_end

```python
def get_intervals_sorted_by_end()
```

Returns list with intervals ordered by their end.

## get\_merged\_intervals

```python
def get_merged_intervals()
```

Returns list of merged intervals so that they do not overlap anymore.

## get\_max\_distance\_between\_intervals

```python
def get_max_distance_between_intervals()
```

Returns the maximum distance between the intervals while accounting for overlap.

## add\_up\_intervals

```python
def add_up_intervals()
```

Returns the added up lengths of intervals without accounting for overlap.

## add\_up\_merged\_intervals

```python
def add_up_merged_intervals()
```

Returns the added up lengths of merged intervals in order to account for overlap.

## get\_thickness

```python
def get_thickness(system, axis)
```

Get the thickness of an atomic system along a basis vector direction.

**Arguments**:

- `system(ase.Atoms)` - The system from which the thickness is evaluated.
- `axis(int)` - The index of the unit cell basis in which direction the
  thickness is evaluated.
  

**Returns**:

- `float` - The thickness of the system as measured from the center of the
  topmost atom to the center of lowermost atom.

## get\_minimized\_cell

```python
def get_minimized_cell(system, axis, min_size)
```

Used to resize the cell in the given cell direction so that all atoms
are within the cell.

**Arguments**:

- `system(ase.Atoms)` - The system to minimize. Must have a cell with
  non-zero volume.
- `axis(int)` - Index of the axis to minimize.
- `axis(min_size)` - The minimum size for the cell in the given direction.
  

**Returns**:

- `ase.Atoms` - The new minimized system.

## cartesian

```python
def cartesian(arrays, out=None)
```

Generate a cartesian product of input arrays.

**Arguments**:

  arrays(sequence of arrays): The arrays from which the product is
  created.
- `out(ndarray)` - Array to place the cartesian product in.
  

**Returns**:

- `ndarray` - 2-D array of shape (M, len(arrays)) containing cartesian
  products formed of input arrays.
  

**Example**:

  --------
  >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
  array([[1, 4, 6],
  [1, 4, 7],
  [1, 5, 6],
  [1, 5, 7],
  [2, 4, 6],
  [2, 4, 7],
  [2, 5, 6],
  [2, 5, 7],
  [3, 4, 6],
  [3, 4, 7],
  [3, 5, 6],
  [3, 5, 7]])

## get\_crystallinity

```python
def get_crystallinity(symmetry_analyser)
```

Quantifies the crystallinity of the structure as a ratio of symmetries
per number of unique atoms in primitive cell. This metric can be used to
distinguish between amorphous and 'regular' crystals.

The number of symmetry operations corresponds to the symmetry operations
corresponding to the hall number of the structure. The symmetry operations
as given by spglib.get_symmetry() are specific to the original structure,
and they have not been reduced to the symmetries of the space group.

**Arguments**:

  symmetry_analyser (): A SymmetryAnalyzer object that has been
  initialized with an atomic structure.
  

**Returns**:

  (float) A ratio of symmetries per unique atoms in the primitive cell.

## get\_distances

```python
def get_distances(system: Atoms, radii="covalent") -> Distances
```

Returns complete distance information.

**Arguments**:

- `system` - The system from which distances are calculated from.

**Returns**:

  A Distances instance.

## get\_radii

```python
def get_radii(radii, atomic_numbers=None) -> np.ndarray
```

Returns an array of atomic radii for each atom.

## swap\_basis

```python
def swap_basis(atoms: Atoms, a: int, b: int)
```

Used to swap two bases in a system. The geometry remains identical.

**Arguments**:

- `a` - First index to swap
- `b` - Second index to swap

## complete\_cell

```python
def complete_cell(a, b, length)
```

Given two basis vectors a and b, creates a third one that is
orthogonal and has the length of 2 * maximum 2D cell height.

**Arguments**:

- `a(np.ndarray)` - First basis vector
- `b(np.ndarray)` - Second basis vector
- `length(float)` - Length of the basis
  

**Returns**:

- `np.ndarray` - The third basis vector
