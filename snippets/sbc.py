import ase.io
from matid.clustering import SBC

system = ase.io.read('system.xyz')

sbc = SBC()
clusters = sbc.get_clusters(system)

for cluster in clusters:
    pass