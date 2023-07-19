from ase.io import read
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.simulation.run import Run, Program
from nomad.datamodel.metainfo.simulation.system import (
    System,
    Atoms as NOMADAtoms
)
from nomad.datamodel.metainfo.simulation.method import (
    Method, BasisSetContainer, BasisSet, Electronic, DFT, XCFunctional, Functional,
    Electronic
)

from nomad.datamodel.metainfo.simulation.calculation import (
    Calculation,
    Energy,
    EnergyEntry,
    BandEnergies
)
from nomad.datamodel.metainfo.simulation.workflow import GeometryOptimization
from nomad.normalizing import normalizers


atoms = read("system.xyz")
template = EntryArchive()
run = template.m_create(Run)
run.program = Program(name='VASP', version='4.6.35')
method = run.m_create(Method)
method.electrons_representation = [BasisSetContainer(
	type='plane waves',
	scope=['wavefunction'],
	basis_set=[BasisSet(
		type='plane waves',
		scope=['valence'],
	)]
)]
method.electronic = Electronic(method='DFT')
xc_functional = XCFunctional(exchange=[Functional(name='GGA_X_PBE')])
method.dft = DFT(xc_functional=xc_functional)
system = run.m_create(System)
system.atoms = NOMADAtoms(
	lattice_vectors=[
		[5.76372622e-10, 0.0, 0.0],
		[0.0, 5.76372622e-10, 0.0],
		[0.0, 0.0, 4.0755698899999997e-10]
	],
	positions=[
		[2.88186311e-10, 0.0, 2.0377849449999999e-10],
		[0.0, 2.88186311e-10, 2.0377849449999999e-10],
		[0.0, 0.0, 0.0],
		[2.88186311e-10, 2.88186311e-10, 0.0],
	],
	labels=['Br', 'K', 'Si', 'Si'],
	periodic=[True, True, True])
scc = run.m_create(Calculation)
scc.system_ref = system
scc.method_ref = method
scc.energy = Energy(
	free=EnergyEntry(value=-1.5936767191492225e-18),
	total=EnergyEntry(value=-1.5935696296699573e-18),
	total_t0=EnergyEntry(value=-3.2126683561907e-22))
template.workflow2 = GeometryOptimization()

template.run[0].calculation[0].system_ref = None
template.run[0].calculation[0].eigenvalues.append(BandEnergies())
template.run[0].calculation[0].eigenvalues[0].kpoints = [[0, 0, 0]]
template.run[0].system = None

system = System()
system.atoms = NOMADAtoms(
	positions=atoms.get_positions() * 1E-10,
	labels=atoms.get_chemical_symbols(),
    atomic_numbers=atoms.get_atomic_numbers(),
	lattice_vectors=atoms.get_cell() * 1E-10,
	periodic=atoms.get_pbc())

template.run[0].m_add_sub_section(Run.system, system)
for normalizer_class in normalizers:
	normalizer = normalizer_class(template)
	normalizer.normalize()

archive = template.m_to_json()
with open("system.archive.json", "w") as fout:
	fout.write(archive)
