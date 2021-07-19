.extra_framework_maker import ExtraFrameworkMaker, ExtraFrameworkAnalyzer
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
from itertools import permutations
from lxml import etree
from contextlib import closing
from collections import OrderedDict


def get_capped_cluster(atoms, file_name):
    """
    Extract smaller cluster containing the extra-framework atoms and then cap all the O.
    Save cluster in both .traj file and .pdb format.
    :param :
    :param :
    :return:
    """
    EFMaker = ExtraFrameworkAnalyzer(atoms)
    cluster = atoms[[index for index in EFMaker.get_extraframework_cluster()]]
    cluster = Zeolite(cluster).cap_atoms()

    write('/Users/jiaweiguo/Box/openMM_test/%s.traj' % file_name, cluster)
    proteindatabank.write_proteindatabank('/Users/jiaweiguo/Box/openMM_test/%s.pdb' % file_name, cluster)
    # gromacs.write_gromacs('/Users/jiaweiguo/Box/cluster.gro', cluster)
    return cluster


def label_pdb(file_name):
    """
    Relabeling the Atom name in proteindatabank file. (required step for openMM)
    The same atom type connecting to different neighboring types are treated differently due to differences in their
    chemical environments, and is therefore named separately.
    """
    filein = open('/Users/jiaweiguo/Box/openMM_test/%s.pdb' % file_name, 'r')
    fileout = open('/Users/jiaweiguo/Box/openMM_test/%s_labeled.pdb' % file_name, 'w')

    name_list = []
    for line in filein.readlines():
        if line.startswith('ATOM') or line.startswith('HETATM'):
            name = line[12:16].strip()
            name_list.append(name)
            name = name + str(name_list.count(name))
            name = name.rjust(4)
            line = line.replace(line[12:16], name, 1)
            # only replacing the first occurrence of line[12:16], atomic symbols are maintained
        fileout.writelines(line)

    filein.close()
    fileout.close()


def get_bonds(cluster, excluded_index=None, excluded_pair=None):
    """
    Using ase.geometry.analysis.Analysis to get all bonds, then remove the repeated ones.
    Function also allows removing certain bonding pair defined by user (excluded_pair).
    Or removing pairs including certain atomic indices (excluded_index).
    :param excluded_index: list of integers
    :param excluded_pair: list of lists
    :return: full bonding list, shortened list.
    If both excluded_index and excluded_pair are None, bonding list == shortened list
    """
    if excluded_index is None:
        excluded_index = []
    if excluded_pair is None:
        excluded_pair = []

    nl = NeighborList(natural_cutoffs(cluster), bothways=True, self_interaction=False)
    nl.update(cluster)

    bond_list, shortened_list = [], []
    for count, indices in enumerate(Analysis(cluster, nl=nl).all_bonds[0]):
        for index in indices:
            if [count, index] not in bond_list and [index, count] not in bond_list:
                bond_list.append([count, index])

    for bond in bond_list:
        if all(single_index not in bond for single_index in excluded_index) and \
                all(tuple(bond) not in list(permutations(pair)) for pair in excluded_pair):
            shortened_list.append(bond)

    return bond_list, shortened_list


def get_angles(cluster, excluded_index=None, excluded_pair=None):
    """
    ase.geometry.analysis.Analysis.unique_angles function does not work, return all angles.
    three-body interactions.
    :param excluded_pair: excluding all [particle1, particle2, particle3] lists involving the excluded pair
    """
    if excluded_index is None:
        excluded_index = []
    if excluded_pair is None:
        excluded_pair = []

    nl = NeighborList(natural_cutoffs(cluster), bothways=True, self_interaction=False)
    nl.update(cluster)

    angle_list, shortened_list = [], []
    for count, indices in enumerate(Analysis(cluster, nl=nl).all_angles[0]):
        for index in indices:
            if all(list(val) not in angle_list for val in list(permutations([count, index[0], index[1]]))):
                angle_list.append([count, index[0], index[1]])

    for angle in angle_list:
        if all(single_index not in angle for single_index in excluded_index) and \
                all(list(value) not in excluded_pair for value in list(permutations(angle, 2))):
            shortened_list.append(angle)

    return angle_list, shortened_list


def get_forces(my_pdb, my_system):
    """
    integrator used for advancing the equations of motion in MD
    doesn't matter what we pick here since we only need the forces on the initial structure, but do need to have it
    :return: forces on atoms in units of eV/A
    """
    integrator = LangevinMiddleIntegrator(3 * kelvin, 1 / picosecond, 0.4 * picoseconds)  # randomly picked
    simulation = Simulation(my_pdb.topology, my_system, integrator)
    simulation.context.setPositions(my_pdb.positions)
    state = simulation.context.getState(getForces=True)
    forces = np.array(state.getForces(asNumpy=True)) * 1.0364e-2 * 0.1  # convert forces from kJ/nm mol to eV/A
    return forces


def get_DFT_forces_single(atoms, atom_index):
    """
    DFT forces for single atoms
    """
    f_vec = atoms.calc.results['forces'][atom_index]  # self.atoms.get_forces()[atom_index]
    f_mag = np.linalg.norm(f_vec)
    return f_vec


def write_xml(atoms, bonds, save_as):
    # on-the-fly generation of force field xml file, matching atoms and bonds with pdb file
    root = etree.Element('ForceField')

    xml_section = etree.SubElement(root, "AtomTypes")
    for atom in atoms:
        properties = {'name': atom.name, 'class': atom.name, 'element': ''.join(filter(lambda x: not x.isdigit(),
                                                                                       atom.name)), 'mass': '0.0'}
        etree.SubElement(xml_section, 'Type', **properties)

    xml_section = etree.SubElement(root, 'Residues')
    xml_residue = etree.SubElement(xml_section, 'Residue', name='MOL')
    for atom in atoms:
        etree.SubElement(xml_residue, 'Atom', name=atom.name, type=atom.name)
    for bond in bonds:
        etree.SubElement(xml_residue, 'Bond', atomName1=bond[0].name, atomName2=bond[1].name)

    tree = etree.ElementTree(root)
    xml = etree.tostring(tree, pretty_print=True).decode('utf-8')

    with closing(open(save_as, 'w')) as f:
        f.write(xml)


def get_FF_forces_single(numb, shortened_bond_list, bond_pair_dict, initial_param_dict, index_to_track=None,
                         shortened_angle_list=None):
    """
    #todo: separate out the force design part
    use numb to keep track of individual configurations
    """
    pdb = PDBFile('/Users/jiaweiguo/Box/openMM_test/cluster_%s_labeled.pdb' % numb)
    atoms = list(pdb.topology.atoms())

    for index in shortened_bond_list:
        pdb.topology.addBond(atoms[index[0]], atoms[index[1]])
    bonds = list(pdb.topology.bonds())

    write_xml(atoms, bonds, '/Users/jiaweiguo/Box/openMM_test/template_test.xml')
    FF = ForceField('/Users/jiaweiguo/Box/openMM_test/template_test.xml')
    system = FF.createSystem(pdb.topology)

    force = add_custom_forces(shortened_bond_list, bond_pair_dict, initial_param_dict)
    system.addForce(force)
    """
    force = HarmonicAngleForce()  # Harmonic angle
    for index in shortened_angle_list:
        force.addAngle(int(index[0]), int(index[1]), int(index[2]), 2.3, 100)
    system.addForce(force)
    """
    if index_to_track is not None:
        return get_forces(pdb, system)[index_to_track]
    else:
        return get_forces(pdb, system)


def add_custom_forces(shortened_bond_list, bond_pair_dict, initial_param_dict):
    force = CustomBondForce("D*(1-exp(-alpha*(r-r0)))^2")  # Morse bond
    force.addPerBondParameter("D")
    force.addPerBondParameter("alpha")
    force.addPerBondParameter("r0")

    for bond in shortened_bond_list:
        for key, ref in bond_pair_dict.items():
            if any(list(val) in ref for val in list(permutations(bond))):
                # print(initial_param_dict[key])
                force.addBond(int(bond[0]), int(bond[1]), initial_param_dict[key])
    return force


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


def func():
    traj = read('/Users/jiaweiguo/Box/MFI_minE_O_less.traj', ':')
    cluster_traj = []
    for count, atoms in enumerate(traj):
        cluster_traj.append(get_capped_cluster(atoms, 'cluster_' + str(count)))
    view(cluster_traj)

    traj = read('/Users/jiaweiguo/Box/MFI_minE_O_less.traj', ':')
    for val in range(len(traj)):
        label_pdb('cluster_%s' % str(val))


def check_atom_types(cluster, index):
    """ assign atom types, same element connected to different neighbors are assigned into different classes.
    For example, extra-framework O (in Cu-O-Cu) is in a different class from framework O (Si-O-Si). Each class
    assignment is unique (each atom belongs to one class and one class only).
    O_EF: extra-framework O
    O-Cu: framework O, connecting to one T-site(Al) and Cu
    O-H: framework O, connecting to one T-site(Al) and H (capping)
    """
    nl = NeighborList(natural_cutoffs(cluster), bothways=True, self_interaction=False)
    nl.update(cluster)

    class_Al = [atom.index for atom in cluster if atom.symbol == 'Al']
    class_Cu = [atom.index for atom in cluster if atom.symbol == 'Cu']
    class_H = [atom.index for atom in cluster if atom.symbol == 'H']
    class_O_EF = [10]
    class_O_Cu = [atom.index for atom in cluster if atom.symbol == 'O' and atom.index not in class_O_EF and
                  all(val not in class_H for val in nl.get_neighbors(atom.index)[0])]
    class_O_H = [atom.index for atom in cluster if atom.symbol == 'O' and atom.index not in class_O_EF + class_O_Cu]
    # print(class_O_Cu)
    # print(class_O_H)

    if index in class_Al:
        return 'class_Al'
    if index in class_Cu:
        return 'class_Cu'
    if index in class_H:
        return 'class_H'
    if index in class_O_EF:
        return 'class_O_EF'
    if index in class_O_Cu:
        return 'class_O_Cu'
    if index in class_O_H:
        return 'class_O_H'
    else:
        return 'None'


def get_bond_or_angle_types(my_list):
    """ assign all bonding pairs or angles into different types based on differences in atom types. For example,
    O(extra-framework)-Cu is different from O(framework)-Cu.
    :param my_list: bond or angle index list of the cluster of interests
    :return type_dict: return a dictionary of all unique bond-pairs or angle types, with "keys" being integers starting
    from 0, and "values" being a list of two atom types string for bonds or three atom types string for angles.
    eg. {0: [AtomClass1, AtomClass2], 1: [AtomClass1, AtomClass3], ...} for bonds
    Note: Bond types such as [AtomClass1, AtomClass2] and [AtomClass2, AtomClass1] are considered the same. Same rules
    also apply for angles.
    :return whole_type_list: return the entire list of bond or angle types assignment of the input.
    len(whole_type_list) = len(my_list)
    """
    type_dict, repeated_list, whole_type_list, count = {}, [], [], 0

    for items in my_list:
        my_list = []
        for val in items:
            my_list.append(check_atom_types(cluster, val))

        whole_type_list.append(my_list)
        if all(list(pair) not in repeated_list for pair in list(permutations(my_list))):
            repeated_list.append(my_list)
            type_dict[count] = my_list
            count += 1
    return type_dict, whole_type_list


def get_index_dict(type_dict, whole_type_list, whole_index_list):
    """ assign bond pairs or angles indices into different bond or angle types, all the pairs or angles within the same
    types will share the same set of force field parameters.
    :param type_dict:
    :param whole_type_list:
    :param whole_index_list:
    :return index_dict: return a dictionary of all bond-pairs or angle indices for each unique bond or angle type,
    using the the same keys as type_dict.
    """
    index_dict = {}
    for key, value in type_dict.items():
        temp_list = []
        for count, items in enumerate(whole_type_list):
            if any(list(pair) == value for pair in list(permutations(items))):
                temp_list.append(whole_index_list[count])
        index_dict[key] = temp_list
    return index_dict


def show_key_value_pair(dict1, dict2):
    """ for better visualization of the bond (or angle) types and bond (or angle) indices that belong to certain types.
    Also displaying the tag (key) of each dict entries, useful for excluding certain types by calling the tags.
    eg. dict1 = bond_type_dict, dict = bond_index_dict
    """
    for key, value in dict1.items():
        print('Tag', key, '-->', value, '-->', dict2[key])


def demo():
    cluster = read('/Users/jiaweiguo/Box/openMM_test/cluster_0.traj', '0')
    bond_list, shortened_bond_list = get_bonds(cluster)
    print(shortened_bond_list == bond_list)

    bond_list, shortened_bond_list = get_bonds(cluster, excluded_index=[2, 3, 8, 9],
                                               excluded_pair=[[11, 12], [0, 4], [0, 6], [1, 5], [1, 7]])
    # for simplicity, excluding all O-H and Al-O bonds
    print(shortened_bond_list)

    angle_list, shortened_angle_list = get_angles(cluster)
    print(shortened_angle_list == angle_list)

    angle_list, shortened_angle_list = get_angles(cluster, excluded_index=[13, 14, 15, 16],
                                                  excluded_pair=[[11, 12], [0, 4], [0, 6], [1, 5], [1, 7]])
    # for simplicity, excluding all angular terms involving O-H and Al-O bonds
    print(shortened_angle_list)

    # bond_list = [[0, 11], [1, 12], [4, 11], [5, 12], [6, 11], [7, 12], [10, 11], [10, 12]]  # for reference

    # load labeled pdb file and build openMM atoms topology
    # pdb file contains only contain basic information such as atom name, symbol, positions, etc
    # can add bonds in pdb file through "CONECT" but unnecessary, can simply use topology.addBond
    numb = 0
    pdb = PDBFile('/Users/jiaweiguo/Box/openMM_test/cluster_%s_labeled.pdb' % numb)
    atoms = list(pdb.topology.atoms())
    print([atom.name for atom in atoms])

    for index in shortened_bond_list:
        pdb.topology.addBond(atoms[index[0]], atoms[index[1]])
    bonds = list(pdb.topology.bonds())
    print([[bond[0].name, bond[1].name] for bond in bonds])
    print([[bond[0].index, bond[1].index] for bond in bonds])

    # force field xml file need a bit attention, but majority of works are one-time-thing (???? maybe not)
    # all pdb files converted from ase.traj are automatically named as MOL for residual name (need to define in the xml
    # file)
    # simplest xml file only need atoms type in residual section, bonds and angles can be added later

    write_xml(atoms, bonds, '/Users/jiaweiguo/Box/openMM_test/template_test.xml')   # on-the-fly generation of ff xml
    FF = ForceField('/Users/jiaweiguo/Box/openMM_test/template_test.xml')
    system = FF.createSystem(pdb.topology)

    # example of customized force function
    force = CustomBondForce("D*(1-exp(-alpha*(r-r0)))^2")  # Morse bond
    force.addPerBondParameter("D")
    force.addPerBondParameter("alpha")
    force.addPerBondParameter("r0")

    # adding bonds in force field, not the same as adding bonds in topology (both are necessary for defining the system)
    for index in shortened_bond_list:
        force.addBond(int(index[0]), int(index[1]), [1.0, 1.0, 2.0])  # need to specify int(), get error otherwise
    # print(force.getNumBonds())
    system.addForce(force)

    # print(get_forces(pdb, system))
    print(get_forces(pdb, system)[10])

    force = HarmonicAngleForce()  # Harmonic angle
    for index in shortened_angle_list:
        force.addAngle(int(index[0]), int(index[1]), int(index[2]), 2.3, 100)
    # print(force.getNumAngles())
    system.addForce(force)

    # print(get_forces(pdb, system))
    print(get_forces(pdb, system)[10])  # predicted forces on O


if __name__ == '__main__':

    """
    # forces on EF-O only
    traj = read('/Users/jiaweiguo/Box/MFI_minE_O_less.traj', ':')
    DFT_f_on_O = []
    for atoms in traj:
        DFT_f_on_O.append(get_DFT_forces_single(traj[0], 288))
    print(DFT_f_on_O)

    FF_f_on_O = []
    for numb in range(len(traj)):
        FF_f_on_O.append(get_FF_forces_single(numb, shortened_bond_list, shortened_angle_list, index_to_track=10))
    print(FF_f_on_O)
    """

    cluster = read('/Users/jiaweiguo/Box/openMM_test/cluster_0.traj', '0')

    bond_list, shortened_bond_list = get_bonds(cluster, excluded_index=[2, 3, 8, 9],
                                               excluded_pair=[[11, 12], [0, 4], [0, 6], [1, 5], [1, 7]])
    # print(bond_list)

    bond_type_dict, whole_bond_type_list = get_bond_or_angle_types(bond_list)
    print(bond_type_dict)
    print(whole_bond_type_list)

    bond_index_dict = get_index_dict(bond_type_dict, whole_bond_type_list, bond_list)
    print(bond_index_dict)

    show_key_value_pair(bond_type_dict, bond_index_dict)
    print('Number of unique bond types:', len(bond_type_dict))

    """
    initial_param_dict = {0: [1, 1, 1], 1: [1.1, 2, 1], 2: [1.2, 4, 2], 3: [1.3, 3, 5], 4: [1.4, 2, 1], 5: [1.5, 1, 1],
                          6: [1.6, 1, 1]}

    # print(get_FF_forces_single(0, shortened_bond_list, bond_index_dict, initial_param_dict))
    """

    """
    angle_list, shortened_angle_list = get_angles(cluster, excluded_index=[13, 14, 15, 16],
                                                  excluded_pair=[[11, 12], [0, 4], [0, 6], [1, 5], [1, 7]])

    angle_type_dict, whole_angle_type_list = get_bond_or_angle_types(angle_list)
    print(angle_type_dict)
    print(whole_angle_type_list)

    angle_index_dict = get_index_dict(angle_type_dict, whole_angle_type_list, angle_list)
    print(angle_index_dict)

    show_key_value_pair(angle_type_dict, angle_index_dict)
    print('Number of unique angle types:', len(angle_type_dict))

    # can be done on shortened list (exclusion based on index)
    angle_type_dict, whole_angle_type_list = get_bond_or_angle_types(shortened_angle_list)
    print(angle_type_dict)
    print(whole_angle_type_list)

    angle_index_dict = get_index_dict(angle_type_dict, whole_angle_type_list, shortened_angle_list)
    print(angle_index_dict)

    show_key_value_pair(angle_type_dict, angle_index_dict)
    print('Number of unique angle types:', len(angle_type_dict))
    """

    # new function, exclude bonds or angles by types
    # starting with bonds

    excluded_type_tag = [0, 1, 3, 6]  # exclude bonds or angles based on the tag number

    shortened_bond_list = []
    for num_key, atom_type_list in bond_type_dict.items():
        if num_key not in excluded_type_tag:
            shortened_bond_list.extend(bond_index_dict[num_key])
    print(shortened_bond_list)
