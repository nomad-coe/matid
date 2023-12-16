---
sidebar_position: 4
---

# Material ID
`SymmetryAnalyzer` can report a 28-digit character sequence for a given structure. This is essentially a hash-digest which is seeded using information about the space group, dimensionality and the Wyckoff sites occupied in the conventional cell. Due to the way MatID generates the conventional cell, this identifier is fairly unique and can be used in many cases to identify structures 3D and 2D materials with "enough" symmetry.

## Usage in search
 The material id can act as a convenient handle for finding materials with very specific structure. As a simple example one could imagine trying to find the forsterite structure. The filters required for this would include:

  - Space group number
  - Wyckoff position occupation (assuming a normalized one is available)

Often this information is not directly reported and if it is, it can be a
daunting task for the user to formulate the queries precisely. This is where the material id can offer a shortcut: one can retrieve it using matid for any structure and the

## Points of failure
This identifier is not guaranteed to be unique, but can be highly useful in practise. In some sense it has a similar value proposal to SMILES for molecules. There are cases in which it will fail:

 - Parametrized Wyckoff sites
 - If isotropic scaling is a problem
 - If variable angles for the conventional cell are a problem
