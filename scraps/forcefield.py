from maze.extra_framework_maker import ExtraFrameworkMaker, ExtraFrameworkAnalyzer
from maze.io_zeolite import read_vasp
from maze.zeolite import PerfectZeolite, Zeolite
from ase.neighborlist import natural_cutoffs, NeighborList
import os
from pathlib import Path
from ase.io import write, read, gromacs, proteindatabank
from ase.visualize import view
import copy
import shutil
from glob import glob
from ase.constraints import FixAtoms
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from ase.geometry.analysis import Analysis
import numpy as np


def get_capped_cluster(atoms, file_name):
    EFMaker = ExtraFrameworkAnalyzer(atoms)
    cluster = atoms[[index for index in EFMaker.get_extraframework_cluster()]]
    cluster = Zeolite(cluster).cap_atoms()

    write('/Users/jiaweiguo/Box/openMM_test/%s.traj' %file_name, cluster)
    proteindatabank.write_proteindatabank('/Users/jiaweiguo/Box/openMM_test/%s.pdb' %file_name, cluster)
    # pdb file labeling needs to be streamlined
    # gromacs.write_gromacs('/Users/jiaweiguo/Box/cluster.gro', cluster)
    return cluster


def get_all_bonds():
    cluster = read('/Users/jiaweiguo/Box/cluster.traj', '0')
    nl = NeighborList(natural_cutoffs(cluster), bothways=True, self_interaction=False)
    nl.update(cluster)

    bond_list, bond_list_wo_H = [], []
    index_H = [a.index for a in cluster if a.symbol == 'H']
    for count, indices in enumerate(Analysis(cluster, nl=nl).all_bonds[0]):
        for index in indices:
            if [count, index] not in bond_list and [index, count] not in bond_list:
                bond_list.append([count, index])
                if count not in index_H and index not in index_H:
                    bond_list_wo_H.append([count, index])

    return bond_list, bond_list_wo_H


def some_random_stuff():
    """
    # example from the OpenMM doc
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    print(simulation.context.getState(getForces=True).getForces(asNumpy=True))
    simulation.reporters.append(PDBReporter('output.pdb', 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
    simulation.step(10000)
    """


if __name__ == '__main__':

    # print(get_all_bonds())

    traj = read('/Users/jiaweiguo/Box/MFI_minE_O_less.traj', ':')
    cluster_traj = []
    for count, atoms in enumerate(traj):
        cluster_traj.append(get_capped_cluster(atoms, 'cluster_'+str(count)))
    view(cluster_traj)








    """
    bond_list = [[0, 11], [1, 12], [4, 11], [5, 12], [6, 11], [7, 12], [10, 11], [10, 12]]
    # shortened, excluding O-H and Al-OH

    pdb = PDBFile('/Users/jiaweiguo/Box/cluster_labeled.pdb')
    atoms = list(pdb.topology.atoms())
    print([atom.name for atom in atoms])

    for index in bond_list:
        pdb.topology.addBond(atoms[index[0]], atoms[index[1]])
    bonds = list(pdb.topology.bonds())
    print([[bond[0].name, bond[1].name] for bond in bonds])

    FF = ForceField('/Users/jiaweiguo/Box/template_v2.xml')   # v2 - no harmonic force functions in xml
    system = FF.createSystem(pdb.topology)

    force = CustomBondForce("D*(1-exp(-alpha*(r-r0)))^2")  # Morse bond
    force.addPerBondParameter("D")
    force.addPerBondParameter("alpha")
    force.addPerBondParameter("r0")
    for index in bond_list:
        force.addBond(index[0], index[1], [1, 1, 2])
    print(force.getNumBonds())
    system.addForce(force)

    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    state = simulation.context.getState(getForces=True)
    forces = np.array(state.getForces(asNumpy=True)) * 1.0364e-2 * 0.1  # convert forces from kJ/nm mol to eV/A
    print(forces)

    angle_list = [11, 10, 12]
    force = HarmonicAngleForce()  # Harmonic angle
    force.addAngle(angle_list[0], angle_list[1], angle_list[2], 2.3, 100)
    print(force.getNumAngles())
    system.addForce(force)

    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    state = simulation.context.getState(getForces=True)
    forces = np.array(state.getForces(asNumpy=True)) * 1.0364e-2 * 0.1  # convert forces from kJ/nm mol to eV/A
    print(forces)
    """
