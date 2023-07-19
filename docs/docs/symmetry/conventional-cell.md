---
sidebar_position: 2
---

# Conventional Cell
For each space group there is a conventional cell that has been standardized by the ITA. spglib can report such cells, but MatID adds an additional step of normalization upon it to make the structure even more standardized.

Although the shape and symmetries of the cell are already well-defined for a conventional cell reported by spglib, the exact location of atoms in the cell can vary. The same structure can have multiple representations, which all look different. If you are only interested in the symmetry properties of the system, this may not be relevant. But if instead you wish to generate a more unique structure for the conventional cell for visualization or identification purposes, then additional processing is required.

The conventional cell reported by MatID makes an additional step of searching for the most unique representation based on chirality-preserving Euclidean Normalizers. These are reporte e.g. by the Bilbao Crystal database. Each normalizer is essentially a transform that can be applied to to the structure to gain a new for of the coventinoal cell without modifying the underlying material. MatID iterates over these different representations and selects one based on the population of different Wyckoff sites which differs among these representations.

It should be noted that certain Wyckoff sites in certain space groups are parametrized. MatID does not perform any selection between these parameters, so in these cases the exact visual representation can still differ. By calling X you can query whether the structure contains parametrized wyckoff positions.

In addition to producing more standardized visualizations for the conventional cell, this step is useful for material identification as seen in the tutorial on Material IDs.