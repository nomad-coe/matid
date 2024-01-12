---
sidebar_position: 1
sidebar_label: SBC Class
---
# SBC Class

Class for performing Symmetry-based clustering (SBC).

You can apply this class for partitioning a larger material into grains, a
heterostructure into it's component etc. The clustering is based on finding
periodically repeating motifs, and as such it is not suitable for e.g.
finding molecules. Any atoms that do not have enough periodic repetitions
will be returned as isolated clusters.

## \_\_init\_\_

```python
def __init__(seed=7)
```

**Arguments**:

- `seed(int)` - The seed that is used for random number generation.

## get\_clusters

```python
def get_clusters(system,
                 angle_tol=20,
                 max_cell_size=6,
                 pos_tol=0.7,
                 merge_threshold=0.5,
                 merge_radius=1,
                 bond_threshold=0.65,
                 overlap_threshold=-0.1,
                 radii="covalent")
```

Used to detect and return structurally separate clusters within the
given system.

**Arguments**:

- `system` _ase.Atoms_ - The structure to partition.
- `angle_tol` _float_ - angle_tol parameter for PeriodicFinder
- `max_cell_size` _float_ - max_cell_size parameter for PeriodicFinder.get_region
- `pos_tol` _float_ - pos_tol parameter for PeriodicFinder.get_region
- `merge_threshold` _float_ - A threshold for merging two clusters
  together. Give as a fraction of shared atoms. Value of 1 would
  mean that clusters are never merged, value of 0 means that they
  are merged always when at least one atom is shared.
- `merge_radius` _float_ - Radius for finding nearby atoms when deciding
  which cluster is closest. The atomic radii are subtracted from
  distances. Given in angstroms.
- `bond_threshold(float)` - Used to control the connectivity threshold
  for defining a chemical connection between atoms. Controls e.g.
  what types of unit cells are accepted and how outliers are
  removed from clusters.
- `overlap_threshold(float)` - Used to exclude non-physical cells by
  checking overlap of atoms. Overlap between two atoms is
  calculated by subtracting atomic radii from the distance between
  the atoms.
- `radii(str|np.ndarray)` - The radii to use for atoms. Use either a preset
  or a custom list of atomic radii for each atom. The available presets are:
  
  - covalent: Covalent radii from DOI:10.1039/B801115J
  - vdw: van Der Waals radii from DOI:10.1039/C3DT50599E
  - vdw_covalent: preferably van Der Waals radii, covalent if vdw
  not defined.
  

**Returns**:

  A list of Clusters.
