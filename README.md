<img src="https://raw.githubusercontent.com/nomad-coe/matid/main/docs/static/img/logo.png" width="300">

![Build status](https://github.com/nomad-coe/matid/actions/workflows/test.yml/badge.svg)
![Coverage Status](./reports/coverage/coverage-badge.svg)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

MatID is a Python package for identifying and analyzing atomistic systems based
on their structure.

# Documentation
For more details and tutorials, visit the documentation at:
[https://nomad-coe.github.io/matid/](https://nomad-coe.github.io/matid/)

You can find even more details in the following articles:

- [Materials structure genealogy and high-throughput topological classification of surfaces and 2D materials](<https://doi.org/10.1038/s41524-018-0107-6>)

## Example: Surface detection and analysis

```python
import ase.io
from ase.visualize import view

from matid.clustering import SBC
from matid.symmetry import SymmetryAnalyzer

# Load structure from a file
system = ase.io.read('data/system.xyz')

# Find interesting substructures using Symmetry-based Clustering (SBC)
sbc = SBC()
clusters = sbc.get_clusters(system)

# Analyze each found cluster printing out the indices of the atoms belonging to
# this cluster and visualizing the conventional cell from which the cluster was
# built from.
for cluster in clusters:

    # Get the indices of the atoms belonging to this cluster
    indices = cluster.indices
    print(indices)

    # Get the dimensionality of the cluster
    dimensionality = cluster.get_dimensionality()
    print(dimensionality)

    # Get the cell from which the cluster is constructed from. The periodicity
    # of this cell indicates in which directions the unit cell has been found to
    # be repeated in (at least once, possibly infinitely).
    cell = cluster.get_cell()
    n_repeated_directions = sum(cell.get_pbc())
    print(n_repeated_directions)

    # Analyze some symmetry properties of the underlying cell to better identify
    # the material from which the cluster has been constructed from.
    analyzer = SymmetryAnalyzer(cell, symmetry_tol=0.5)
    conv_sys = analyzer.get_conventional_system()
    view(conv_sys)
```

# Installation

## pip
```sh
pip install matid
```

## From source
```sh
git clone https://github.com/nomad-coe/matid.git
cd matid
pip install .
```

