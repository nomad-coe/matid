from matid.clustering import SBC
import ase.io
from ase.visualize import view

system = ase.io.read('system.xyz')

sbc = SBC()
clusters = sbc.get_clusters(system)

for cluster in clusters:
    print(cluster.indices)
    view(cluster.get_cell())
