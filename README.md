<img src="https://raw.githubusercontent.com/nomad-coe/matid/develop/docs/static/img/logo.png" width="300">

![Build status](https://github.com/nomad-coe/matid/actions/workflows/test.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/nomad-coe/matid/badge.svg?branch=main)](https://coveralls.io/github/nomad-coe/matid?branch=main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

MatID is a Python package for identifying and analyzing atomistic systems based
on their structure.

# Documentation
For more details and tutorials, visit the documentation at:
[https://nomad-coe.github.io/matid/](https://nomad-coe.github.io/matid/)

You can find even more details in the following articles:

- [Materials structure genealogy and high-throughput topological classification of surfaces and 2D materials](<https://doi.org/10.1038/s41524-018-0107-6>)

## Example: Surface detection and analysis

```python
import numpy as np
from ase.visualize import view
from ase.build import bcc100, molecule
from matid import Classifier, SymmetryAnalyzer

# Generating a surface adsorption geometry with ASE.
adsorbent = bcc100('Fe', size=(3, 3, 4), vacuum=8)
adsorbate = molecule("H2O")
adsorbate.rotate(180, [1, 0, 0])
adsorbate.translate([4.3, 4.3, 13.5])
system = adsorbent + adsorbate
system.set_pbc([True, True, True])

# Add noise and defects to the structure
positions = system.get_positions()
positions += 0.25*np.random.rand(*positions.shape)
system.set_positions(positions)
del system[31]

# Visualize the final system
view(system)

# Run the classification
classifier = Classifier(pos_tol=1.0, max_cell_size=6)
classification = classifier.classify(system)

# Print classification
print("Structure classified as: {}".format(classification))

# Print found outliers
outliers = classification.outliers
print("Outlier atoms indices: {}".format(outliers))

# Visualize the cell that was found by matid
prototype_cell = classification.prototype_cell
view(prototype_cell)

# Visualize the corresponding conventional cell
analyzer = SymmetryAnalyzer(prototype_cell, symmetry_tol=0.5)
conv_sys = analyzer.get_conventional_system()
view(conv_sys)

# Visualize the corresponding primitive cell
prim_sys = analyzer.get_primitive_system()
view(prim_sys)

# Print space group number
spg_number = analyzer.get_space_group_number()
print("Space group number: {}".format(spg_number))
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

