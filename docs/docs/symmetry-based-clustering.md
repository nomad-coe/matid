---
sidebar_position: 2
---

# Symmetry-based Clustering (SBC)

There are several tools for building atomistic systems from components, but the
reverse process of identifying these components from an existing system is much
harder. MatID utilizes a custom clustering algorithm, called Symmetry-based
Clustering (SBC), which can cluster atoms based on local translational symmetry.

Clustering is performed using the `SBC` class. The basic syntax is relatively
simple: you have to initialize the SBC class with parameters that are suitable
for your use case (sensible defaults are provided), and then call the
`get_clusters()`-method. The following demonstrates this on a existing structure
file:

```python
from ase.io import read
from matid.clustering import Clusterer

system = ase.io.read('system.xyz')

sbc = SBC()
clusters = SBC.get_clusters(system)
```

The return value is a list of the found `Cluster` instances. You can perform
further analysis on these found clusters by using methods and attributes of this
class. E.g. to visualize all of the found clusters, we could do the following:

```python
for cluster in clusters:

```

The resulting image is as follows: