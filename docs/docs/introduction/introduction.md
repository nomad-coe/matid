---
sidebar_position: 1
---

# Introduction
MatID is a Python package for identifying and analyzing atomistic systems based on their structure. MatID is designed to help researchers in the automated analysis and labeling of atomistic data.

The package comes with different tools that can be connected together to form a
workflow for analyzing structures. Here is a brief list of the different tools
and their functionality:

- Clustering: Used to identify and report repeating structures in the given atomistic system. Works both with finite and periodic systems.
- Dimensionality analysis: Used to determine the intended dimensionality of a system taking periodic boundary conditions into account.
- Symmetry analysis: Routines for analyzing the symmetry of structures. Can be used to return a highly unique conventional cell, detailed Wyckoff position information etc.

It is up to the user to decide how to connect these tools in order to analyze a given system. For example, If one wanted to identify components for a heterostructure containing layered materials, we could do something like this:

```python
# Load system to analyze

# Get the dimensionality of the system

# Perform clustering

# Retrieve conventional cell and a material identifier for each detected cluster

```