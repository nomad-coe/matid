---
sidebar_position: 4
---

# Dimensionality Analysis
MatID can be used to find out the dimensionality of an atomistic geometry. The
definition of dimensionality can be quite tricky, so it will be good to know the
basics behind our definition. Dimensionality in MatID could be defined as: __the number of dimensions in which the structre is connected to itself through periodic boundary conditions__. What is meant by **connected** is explained below in section X. It should be noted that MatID does **not consider the shape or size** of a structure when assigning dimensionality, but instead uses the periodic boundary conditions of your system to determine the intended dimensionality.

To determine the dimensionality of a system, MatID uses a modified version of the
topological scaling algorithm (TSA) [1]_. The algorithm is based on analyzing
the size scaling of atomic clusters when going from the original system to a
bigger supercell of the same system. With TSA, the dimensionality $D$ is
given by

$$
	D=\begin{cases}
		n_\text{pbc}-\log_n (N_{n}) \text{, when}~n_\text{pbc} \neq 0  \\
		0\text{, when}~n_\text{pbc} = 0
	\end{cases}
$$

where $N_n$ is the number of clusters in a supercell that is repeated
$n$ times in each periodic direction and $n_\mathrm{pbc}$ is the
number of periodic dimensions. For the clustering we use the Density-Based
Spatial Clustering of Applications with Noise [2]_ data clustering algorithm.
The advantage of this algorithm is that it does not require an initial guess
for the number of clusters and it can find arbitrarily shaped clusters. The
clustering requires that we define a metric for the distance between atoms. We
use the following metric:

$$
   d_{ij} = \lvert \vec{r}_i - \vec{r}_j \rvert^{\text{MIC}} - r_i - r_j
$$

where $\vec{r}_i$ and $\vec{r}_i$ are the cartesian positions of
atom $i$ and $j$, respectively, and $r_i$ and $r_j$ are
their radii. The radii definition can be changed and defaults to covalent radii
[3]_ . It is important to notice that in this metric the distances always
follow the minimum image convention (MIC), i.e.  the distance is calculated
between two closest periodic neighbours. By using the distance to the closest
periodic neighbour we obtain the correct clusters regardless of what shape of
cell is used in the original simulation.

The clustering uses two parameters: the minimum cluster size
$n_\mathrm{min}$ and the neighbourhood radius $\epsilon$. We set
$n_\mathrm{min}$ to 1 to allow clusters consisting of even single atoms
and $\epsilon$ defaults to 3.5 Å. At present, a system, in which there is
more than one cluster in the original non-repeated system ($N_1 \gt 1$),
is classified as unknown. Such a case corresponds to systems with multiple
components that are spatially separated, such as a molecule far above a
surface, low density gases, widely spaced clusters in vacuum, etc.

The following code illustrates how dimensionality detection can be performed
with MatID.

.. literalinclude:: ../../../examples/dimensionality.py
   :language: python

This example if also available in "examples/dimensionality.py".

.. [1] Ashton, M., Paul, J., Sinnott, S. B. & Hennig, R. G. Topology-scaling identification of layered solids and stable exfoliated 2d materials. Phys. Rev. Lett. 118, 106101 (2017)

.. [2] Ester, M., Kriegel, H.-P., Sander, J. & Xu, X. A density-based algorithm for discovering clusters in large spatial databases with noise. KDD’96 Proceedings of the Second International Conference on Knowledge Discovery and Data Mining 226–231 (1996).

.. [3] Cordero, B. et al. Covalent radii revisited. Dalton Trans. 2832–2838 (2008)