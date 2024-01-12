---
sidebar_position: 2
sidebar_label: Cluster Class
---
# Cluster Class

Contains information about a cluster, i.e. a part of some bigger original
atomistic system.

This class is a simple data container where every attribute should be set
only once, but does not have to be set in one go and can be insted built
gradually.

## \_\_init\_\_

```python
def __init__(indices=None,
             species=None,
             region=None,
             dimensionality=None,
             cell=None,
             system=None,
             distances=None,
             radii=None,
             bond_threshold=None)
```

**Arguments**:

- `indices(Iterable)` - Contains the indices of atoms belonging to this
  cluster.
- `species(set)` - Contains the species of atoms belonging to this
  cluster. Each unique species should be include only once.
- `region(Region)` - The Region instance from which this cluster was
  exracted from.
- `dimensionality(int)` - The dimensionality of the cluster. Can be set
  initially here or calculated through the get_dimensionality-function.
- `cell(ase.Atoms)` - The unit cell from which this cluster is
  constructed from.
- `system(ase.Atoms)` - Reference to the original system which this
  cluster is a part of.
- `distances(Distances)` - Contains cached distance information about
  this cluster.
- `radii(ndarray)` - Contains the radii for each atom in the cluster as
  floating point numbers.

## get\_cell

```python
def get_cell() -> Atoms
```

Used to fetch the prototypical cell for this cluster if one exists.

## get\_atoms

```python
def get_atoms() -> Atoms
```

Returns the ase.Atoms object for this cluster.

## get\_dimensionality

```python
def get_dimensionality() -> int
```

Shortcut for fetching the dimensionality of the cluster using
matid.geometry.get_dimensionality and the radii + bond thresholds that
were used during the clustering.
