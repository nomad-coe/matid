import numpy as np
import ase.io
from ase.visualize import view

from matid.clustering import SBC
from matid.symmetry import SymmetryAnalyzer
from matid.geometry import get_dimensionality

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
	dimensionality = get_dimensionality(system[indices])
	print(dimensionality)

	# Visualize the conventional cell which the cluster is based on
	analyzer = SymmetryAnalyzer(cluster.cell(), symmetry_tol=0.5)
	conv_sys = analyzer.get_conventional_system()
	view(conv_sys)
