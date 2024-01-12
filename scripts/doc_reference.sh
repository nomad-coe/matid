# Used to build the documentation using PyDoc:
# https://niklasrosenstein.github.io/pydoc-markdown/

# SBC Class
cd ..
pydoc-markdown -m matid.clustering.sbc docs/pydoc-markdown.yml > docs/docs/reference/sbc.md
echo "---
sidebar_position: 1
sidebar_label: SBC Class
---
$(cat docs/docs/reference/sbc.md)" > docs/docs/reference/sbc.md

# Cluster Class
pydoc-markdown -m matid.clustering.cluster docs/pydoc-markdown.yml > docs/docs/reference/cluster.md
echo "---
sidebar_position: 2
sidebar_label: Cluster Class
---
$(cat docs/docs/reference/cluster.md)" > docs/docs/reference/cluster.md

# SymmetryAnalyzer Class
pydoc-markdown -m matid.symmetry.symmetryanalyzer docs/pydoc-markdown.yml > docs/docs/reference/symmetryanalyzer.md
echo "---
sidebar_position: 3
sidebar_label: SymmetryAnalyzer Class
---
$(cat docs/docs/reference/symmetryanalyzer.md)" > docs/docs/reference/symmetryanalyzer.md

# Geometry Module
pydoc-markdown -m matid.geometry.geometry docs/pydoc-markdown.yml > docs/docs/reference/geometry.md
echo "---
sidebar_position: 4
sidebar_label: Geometry Module
---
# Geometry Module
$(cat docs/docs/reference/geometry.md)" > docs/docs/reference/geometry.md