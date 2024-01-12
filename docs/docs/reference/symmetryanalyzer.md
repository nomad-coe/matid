---
sidebar_position: 3
sidebar_label: SymmetryAnalyzer Class
---
# SymmetryAnalyzer Class

A base class for getting symmetry related properties of unit cells.

## \_\_init\_\_

```python
def __init__(system=None, symmetry_tol=None, min_2d_thickness=1)
```

**Arguments**:

- `system(ASE.Atoms)` - The system to inspect.
- `symmetry_tol(float)` - The tolerance for the symmetry detection.
- `min_2d_thickness(float)` - The minimum thickness in angstroms for the
  conventional cell that is returned for 2D systems.

## set\_system

```python
def set_system(system)
```

Sets a new system for analysis.

## reset

```python
def reset()
```

Used to reset all the cached values.

## get\_material\_id

```python
def get_material_id()
```

Returns a 28-character identifier for this material. The identifier
is calculated by hashing a set of the symmetry properties found in the
material, including:

- Space group number
- Wyckoff position letters and the species occupied in them

## get\_space\_group\_number

```python
def get_space_group_number()
```

**Returns**:

- `int` - The space group number.

## get\_space\_group\_international\_short

```python
def get_space_group_international_short()
```

**Returns**:

- `str` - The international space group short symbol.

## get\_hall\_symbol

```python
def get_hall_symbol()
```

**Returns**:

- `str` - The Hall symbol.

## get\_hall\_number

```python
def get_hall_number()
```

**Returns**:

- `int` - The Hall number.

## get\_point\_group

```python
def get_point_group()
```

Symbol of the crystallographic point group in the Hermann-Mauguin
notation.

**Returns**:

- `str` - point group symbol

## get\_is\_chiral

```python
def get_is_chiral()
```

Returns a boolean value that tells if this object is chiral or not
(achiral). A chiral object has symmetry operations that are all proper,
i.e. their determinant is +1.

**Returns**:

- `bool` - is the object chiral.

## get\_has\_free\_wyckoff\_parameters

```python
def get_has_free_wyckoff_parameters()
```

Tells whether this system has Wyckoff positions with free variables.

**Returns**:

- `bool` - Indicates the presence of Wyckoff positions with free variables.

## get\_crystal\_system

```python
def get_crystal_system()
```

Get the crystal system based on the space group number. There are
seven different crystal systems:

- Triclinic
- Monoclinic
- Orthorhombic
- Tetragonal
- Trigonal
- Hexagonal
- Cubic

**Returns**:

- `str` - The name of the crystal system.

## get\_bravais\_lattice

```python
def get_bravais_lattice()
```

Return Bravais lattice in the Pearson notation, where the first
lowercase letter indicates the crystal system, and the second uppercase
letter indicates the centring type.

Crystal system letters:

- a = triclinic
- m = monoclinic
- o = orthorhombic
- t = tetragonal
- h = hexagonal and trigonal
- c = cubic

Lattice type letters:

- P = Primitive
- S (= A or B or C) = One side/face centred
- I = Body centered
- R = Rhombohedral centring
- F = All faces centred

:param crystal_system: The crystal system
:param space_group: The space group number.
:type crystal_system: str
:type space_group: int

:return: The Bravais lattice in the Pearson notation.
:rtype: str

## get\_primitive\_system

```python
def get_primitive_system()
```

Returns a primitive description for this system.

This description uses a primitive lattice where positions of the
atoms, and the cell basis vectors are idealized to follow the
symmetries that were found with the given precision. This means that
e.g. the volume, density, angles between basis vectors and basis vector
lengths may have small deviations from the original system.

**Returns**:

- `ASE.Atoms` - The primitive system.

## get\_conventional\_system

```python
def get_conventional_system()
```

Used to get the conventional representation of this system.

This description uses a conventional lattice where positions of the
atoms, and the cell basis vectors are idealized to follow the
symmetries that were found with the given precision. This means that
e.g. the volume, density, angles between basis vectors and basis vector
lengths may have small deviations from the original system.

## get\_rotations

```python
def get_rotations()
```

Get the rotational parts of the Seitz matrices that are associated
with this space group. Each rotational matrix is accompanied by a
translation with the same index.

**Returns**:

- `np.ndarray` - Rotation matrices.

## get\_translations

```python
def get_translations()
```

Get the translational parts of the Seitz matrices that are
associated with this space group. Each translation is accompanied
by a rotational matrix with the same index.

**Returns**:

- `np.ndarray` - Translation vectors.

## get\_choice

```python
def get_choice()
```

**Returns**:

- `str` - A string specifying the centring, origin and basis vector
  settings.

## get\_wyckoff\_letters\_original

```python
def get_wyckoff_letters_original()
```

**Returns**:

  list of str: Wyckoff letters for the atoms in the original system.

## get\_equivalent\_atoms\_original

```python
def get_equivalent_atoms_original()
```

The equivalent atoms are the same as what spglib already outputs, as
changes in the wyckoff letters will not afect the equivalence:

**Returns**:

  list of int: A list that maps each atom into a symmetry equivalent
  set.

## get\_wyckoff\_letters\_conventional

```python
def get_wyckoff_letters_conventional()
```

Get the Wyckoff letters of the atoms in the conventional system.

**Returns**:

  list of str: Wyckoff letters.

## get\_wyckoff\_sets\_conventional

```python
def get_wyckoff_sets_conventional(return_parameters=True)
```

Get a list of Wyckoff sets for this system. Wyckoff sets combine
information about the atoms and their positions at specific Wyckoff
positions.

**Arguments**:

- `return_parameters` _bool_ - Whether to return the value of possible
  free Wyckoff parameters. Set to false if they are not needed,
  as their determination can take some time.
  

**Returns**:

  list of WyckoffSets: A list of :class:`.WyckoffSet` objects for the
  conventional system.

## get\_equivalent\_atoms\_conventional

```python
def get_equivalent_atoms_conventional()
```

List of equivalent atoms in the idealized system.

**Returns**:

  list of int: A list that maps each atom into a symmetry equivalent
  set.

## get\_wyckoff\_letters\_primitive

```python
def get_wyckoff_letters_primitive()
```

Get the Wyckoff letters of the atoms in the primitive system.

**Returns**:

  list of str: Wyckoff letters.

## get\_equivalent\_atoms\_primitive

```python
def get_equivalent_atoms_primitive()
```

List of equivalent atoms in the primitive system.

**Returns**:

  list of int: A list that maps each atom into a symmetry equivalent
  set.

## get\_symmetry\_dataset

```python
def get_symmetry_dataset()
```

Calculates the symmetry dataset with spglib for the given system.

## get\_symmetry\_operations

```python
def get_symmetry_operations()
```

The symmetry operations of the original structure as rotations and
translations.

**Returns**:

  Dictionary containing an entry for rotations containing a np.array
  with 3x3 matrices for each symmetry operation and an entry
  "translations" containing np.array of translations for each
  symmetry operation.
  3*1 np.ndarray: The shift of the origin as a vector.
