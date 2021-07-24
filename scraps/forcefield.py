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
from itertools import permutations
from lxml import etree
from contextlib import closing
from collections import OrderedDict
from scipy.optimize import least_squares, minimize
import matplotlib.pyplot as plt


def get_capped_cluster(atoms, file_name):
    """ #TODO: check whether capping is necessary
    Extract smaller cluster containing the extra-framework atoms and then cap all the O.
    Save cluster in both .traj file and .pdb format.
    :param :
    :param :
    :return:
    """
    EFMaker = ExtraFrameworkAnalyzer(atoms)
    cluster = atoms[[index for index in EFMaker.get_extraframework_cluster()]]
    # cluster = Zeolite(cluster).cap_atoms()

    write('/Users/jiaweiguo/Box/openMM_test/%s.traj' % file_name, cluster)
    proteindatabank.write_proteindatabank('/Users/jiaweiguo/Box/openMM_test/%s.pdb' % file_name, cluster)
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


def get_bonds(cluster, mult=1, excluded_index=None, excluded_pair=None):
    """
    Using ase.geometry.analysis.Analysis to get all bonds, then remove the repeated ones.
    Function also allows removing certain bonding pair defined by user (excluded_pair).
    Or removing pairs including certain atomic indices (excluded_index).
    :param cluster:
    :param mult:
    :param excluded_index: list of integers
    :param excluded_pair: list of lists
    :return: full bonding list, shortened list.
    If both excluded_index and excluded_pair are None, bonding list == shortened list
    """
    if excluded_index is None:
        excluded_index = []
    if excluded_pair is None:
        excluded_pair = []

    nl = NeighborList(natural_cutoffs(cluster, mult=mult), bothways=True, self_interaction=False)
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


def get_angles(cluster, mult=1, excluded_index=None, excluded_pair=None):
    """
    #TODO: consider combining get_bonds and get_angles function
    ase.geometry.analysis.Analysis.unique_angles function does not work, return all angles.
    three-body interactions.
    :param excluded_pair: excluding all [particle1, particle2, particle3] lists involving the excluded pair
    """
    if excluded_index is None:
        excluded_index = []
    if excluded_pair is None:
        excluded_pair = []

    nl = NeighborList(natural_cutoffs(cluster, mult=mult), bothways=True, self_interaction=False)
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
        return 'Al'
    if index in class_Cu:
        return 'Cu'
    if index in class_H:
        return 'H'
    if index in class_O_EF:
        return 'O-EF'
    if index in class_O_Cu:
        return 'O-Cu'
    if index in class_O_H:
        return 'O-H'
    else:
        return 'None'


def get_property_types(cluster, property_list):
    """ assign all bonding pairs or angles into different types based on differences in atom types. For example,
    O(extra-framework)-Cu is different from O(framework)-Cu.
    :param property_list: bond or angle index list of the cluster of interests
    :return type_dict: return a dictionary of all unique bond-pairs or angle types, with "keys" being integers starting
    from 0, and "values" being a list of two atom types string for bonds or three atom types string for angles.
    eg. {0: [AtomClass1, AtomClass2], 1: [AtomClass1, AtomClass3], ...} for bonds
    Note: Bond types such as [AtomClass1, AtomClass2] and [AtomClass2, AtomClass1] are considered the same. Same rules
    also apply for angles.
    :return whole_type_list: return the entire list of bond or angle types assignment of the input.
    len(whole_type_list) = len(my_list)
    """
    type_dict, repeated_list, whole_type_list, count = {}, [], [], 0

    for items in property_list:
        my_list = []
        for val in items:
            my_list.append(check_atom_types(cluster, val))
        whole_type_list.append(my_list)
        if all(list(pair) not in repeated_list for pair in list(permutations(my_list))):
            repeated_list.append(my_list)
            type_dict[count] = my_list
            count += 1

    return type_dict, whole_type_list


def _get_index_dict(type_dict, whole_type_list, index_list):
    """ assign bond pairs or angles indices into different bond or angle types, all the pairs or angles within the same
    types will share the same set of force field parameters.
    :param type_dict:
    :param whole_type_list:
    :param index_list:
    :return index_dict: return a dictionary of all bond-pairs or angle indices for each unique bond or angle type,
    using the the same keys as type_dict.
    """
    index_dict = {}
    for key, value in type_dict.items():
        temp_list = []
        for count, items in enumerate(whole_type_list):
            if any(list(pair) == value for pair in list(permutations(items))):
                temp_list.append(index_list[count])
        index_dict[key] = temp_list

    return index_dict


def get_type_index_pair(type_dict, whole_type_list, index_list):
    """ write bond_type and bond_index into a single dictionary; can use tuples as dictionary key, not lists
    :param type_dict:
    :param whole_type_list:
    :param index_list:
    """
    bond_index_dict = _get_index_dict(type_dict, whole_type_list, index_list)
    type_index_dict = {}
    for key, value in type_dict.items():
        type_index_dict[tuple(value)] = bond_index_dict[key]
    return type_index_dict


def pretty_print(my_dict):
    """ for better visualization of the bond (or angle) types and bond (or angle) indices that belong to certain types.
    """
    for key, value in my_dict.items():
        print(key, '-->', value)


def shorten_index_list_by_types(type_index_dict, exclude_atom_type=None, exclude_property_type=None,
                                include_property_type=None, case=0):
    """
    allow excluding certain property types or only including certain types
    """

    if exclude_atom_type is not None and exclude_property_type is None:
        case = 1
    if exclude_property_type is not None and exclude_atom_type is None:
        case = 2
    if exclude_property_type is not None and exclude_atom_type is not None:
        case = 3
    if include_property_type is not None:
        case = 4

    shortened_list = []
    for type_list, index_list in type_index_dict.items():
        if case == 1 and all(single_type not in type_list for single_type in exclude_atom_type):
            shortened_list.extend(index_list)
        elif case == 2 and all(list(value) not in exclude_property_type for value in list(permutations(type_list))):
            shortened_list.extend(index_list)
        elif case == 3 and all(single_type not in type_list for single_type in exclude_atom_type) and \
                all(list(value) not in exclude_property_type for value in list(permutations(type_list))):
            shortened_list.extend(index_list)
        elif case == 4 and any(list(value) in include_property_type for value in list(permutations(type_list))):
            shortened_list.extend(index_list)

    return shortened_list


def set_up_openMM_system(numb, shortened_bond_list):
    """ Feed pdb topology file and xml force field file into openMM, generate a system for the MD simulation/force
    calculation.
    :return pdb:
    :return system:
    """
    pdb = PDBFile('/Users/jiaweiguo/Box/openMM_test/cluster_%s_labeled.pdb' % numb)
    atoms = list(pdb.topology.atoms())

    for index in shortened_bond_list:
        pdb.topology.addBond(atoms[index[0]], atoms[index[1]])
    bonds = list(pdb.topology.bonds())

    write_xml(atoms, bonds, '/Users/jiaweiguo/Box/openMM_test/template_test.xml')
    FF = ForceField('/Users/jiaweiguo/Box/openMM_test/template_test.xml')
    system = FF.createSystem(pdb.topology)
    return pdb, system


def custom_openMM_force_object(system, bond_list, bond_type_index_dict, bond_param_dict, angle_list=None,
                               angle_type_index_dict=None, angle_param_dict=None):
    """ #todo: make this function for flexible for addition/manipulation
    :param bond_list: list to be included into force field
    :param angle_list:
    :param bond_type_index_dict: {(type): [index], ...}
    :param angle_type_index_dict:
    :param bond_param_dict: {(type): [param], ...} Note: parameters here uses the standard units, kJ, nm, ...
    :param angle_param_dict:
    :return system: openMM system with custom forces added onto it
    """
    force = CustomBondForce("D*(1-exp(-alpha*(r-r0)))^2")  # Morse bond
    force.addPerBondParameter("D")
    force.addPerBondParameter("alpha")
    force.addPerBondParameter("r0")

    for bond in bond_list:
        for my_type, my_index in bond_type_index_dict.items():
            if any(list(val) in my_index for val in list(permutations(bond))):
                force.addBond(int(bond[0]), int(bond[1]), bond_param_dict.get(my_type))
    system.addForce(force)

    force = HarmonicAngleForce()  # Harmonic angle
    for angle in angle_list:
        for my_type, my_index in angle_type_index_dict.items():
            if any(list(val) in my_index for val in list(permutations(angle))):
                type_tag = [tuple(val) for val in list(angle_param_dict.keys()) if val in list(permutations(my_type))]
                force.addAngle(int(angle[0]), int(angle[1]), int(angle[2]), *angle_param_dict.get(type_tag[0]))
    system.addForce(force)

    return system


def get_FF_forces(pdb, system, bond_list, bond_type_index_dict, bond_param_dict, angle_list=None,
                               angle_type_index_dict=None, angle_param_dict=None):
    """
    use numb to keep track of individual configurations
    integrator used for advancing the equations of motion in MD
    doesn't matter what we pick here since we only need the forces on the initial structure, but do need to have it
    :return: forces values on atoms in units of eV/A
    """
    system = custom_openMM_force_object(system, bond_list, bond_type_index_dict, bond_param_dict, angle_list,
                               angle_type_index_dict, angle_param_dict)

    integrator = LangevinMiddleIntegrator(3 * kelvin, 1 / picosecond, 0.4 * picoseconds)  # randomly picked
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    state = simulation.context.getState(getForces=True)
    forces = np.array(state.getForces(asNumpy=True)) * 1.0364e-2 * 0.1  # convert forces from kJ/nm mol to eV/A

    return forces


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


# section below deals with multiple input structures for force field training
def prep_topologies():
    """
    Save clusters into .traj as well, for later comparison and trouble shooting
    """
    traj = read('/Users/jiaweiguo/Box/MFI_minE_O_less.traj', ':')
    cluster_traj = []
    for count, atoms in enumerate(traj):
        cluster_traj.append(get_capped_cluster(atoms, 'cluster_' + str(count)))
    # view(cluster_traj)

    for val in range(len(traj)):
        label_pdb('cluster_%s' % str(val))


def get_DFT_forces_single(atoms, atom_index):
    """
    reference DFT forces on single atoms
    """
    f_vec = atoms.calc.results['forces'][atom_index]  # self.atoms.get_forces()[atom_index]
    f_mag = np.linalg.norm(f_vec)
    return f_vec


def get_required_objects_for_ff(numb, included_bond_type, included_angle_type):
    """ To reduce computational cost, objects such as pdb, system, shortened_bond_list, bond_type_index_dict are kept
    fixed for each configuration during the optimization (only run once).
    """
    cluster = read('/Users/jiaweiguo/Box/openMM_test/cluster_%s.traj' % numb, '0')

    bond_index_list, shortened_bond_index_list = get_bonds(cluster, mult=2)
    bond_type_dict, whole_bond_type_list = get_property_types(cluster, bond_index_list)
    # print('Number of unique bond types:', len(bond_type_dict))

    angle_index_list, shortened_angle_index_list = get_angles(cluster, mult=2)
    angle_type_dict, whole_angle_type_list = get_property_types(cluster, angle_index_list)
    # print('Number of unique angle types:', len(angle_type_dict))

    bond_type_index_dict = get_type_index_pair(bond_type_dict, whole_bond_type_list, bond_index_list)
    # pretty_print(bond_type_index_dict)

    angle_type_index_dict = get_type_index_pair(angle_type_dict, whole_angle_type_list, angle_index_list)
    # pretty_print(angle_type_index_dict)

    shortened_bond_list = shorten_index_list_by_types(bond_type_index_dict, include_property_type=included_bond_type)
    shortened_angle_list = shorten_index_list_by_types(angle_type_index_dict, include_property_type=included_angle_type)
    # print(shortened_bond_list)

    pdb, system = set_up_openMM_system(numb, shortened_bond_list)

    return pdb, system, shortened_bond_list, bond_type_index_dict, shortened_angle_list, angle_type_index_dict


def make_parity_plot(ff_forces, dft_forces):
    plt.figure()
    fig, ax = plt.subplots()
    plt.plot(dft_forces, ff_forces, 'o')
    plt.xlabel('DFT_force', fontsize=18)
    plt.ylabel('FF_force', fontsize=18)
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    plt.title('Force fitting on O', fontsize=18)
    plt.show()


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

    # save openMM properties of each configuration into dict
    traj = read('/Users/jiaweiguo/Box/MFI_minE_O_less.traj', ':')
    included_bond_type, included_angle_type = [['Al', 'Cu'], ['O-Cu', 'Cu'], ['O-EF', 'Cu']], [['Cu', 'O-EF', 'Cu']]
    #included_bond_type, included_angle_type = [['O-EF', 'Cu']], [['Cu', 'O-EF', 'Cu']]

    info_dict = {}
    for numb in range(len(traj)):
        pdb, system, shortened_bond_list, bond_type_index_dict, shortened_angle_list, angle_type_index_dict = \
            get_required_objects_for_ff(numb, included_bond_type, included_angle_type)
        info_dict[numb] = [pdb, system, shortened_bond_list, bond_type_index_dict, shortened_angle_list,
                           angle_type_index_dict]
    # print(info_dict)

    param1 = [ 3.20096548e+02, 5.29980417e+00, 8.27081142e-02, -8.23960437e+01, 2.05108589e+00, 2.78910606e-01,
               -4.79268191e+02, 3.99067733e+00, 2.01045083e-01]
    # [1.2, 4, 0.2, 1.2, 4, 0.2, 1.2, 4, 0.2]
    param2 = [1.51297586e-01, 3.53318807e+01]  # [2.3, 40]
    # get force field forces
    bond_param_dict = {('Al', 'Cu'): param1[0:3], ('O-Cu', 'Cu'): param1[0:3], ('O-EF', 'Cu'): param1[0:3]}
    #bond_param_dict = {('O-EF', 'Cu'): param1}
    angle_param_dict = {('Cu', 'O-EF', 'Cu'): param2}  # angle:radians, k:kJ/mol/radian^2

    ff_f_on_O = []
    my_dict = copy.deepcopy(info_dict)
    for config_tag, info_list in my_dict.items():
        ff_forces = get_FF_forces(info_list[0], info_list[1], info_list[2], info_list[3], bond_param_dict, info_list[4],
                                  info_list[5], angle_param_dict)[10:13]
        ff_f_on_O.append([force_list for force_list in ff_forces])
    print(ff_f_on_O)

    DFT_f_on_O = []
    for atoms in traj:
        DFT_f_on_O.append([get_DFT_forces_single(atoms, atom_index=val) for val in [288, 289, 290]])
    print(DFT_f_on_O)

    make_parity_plot(np.array(np.reshape(ff_f_on_O, [-1, 3])), np.array(np.reshape(DFT_f_on_O, [-1, 3])))

    """
    def get_residue(param, info_dict):
        # todo: rewrite into class, self.included_bond_type and self.info_dict are better

        bond_param_dict = {tuple(included_bond_type[0]): list(param[0:3]),
                           tuple(included_bond_type[1]): list(param[3:6]),
                           tuple(included_bond_type[2]): list(param[6:9])}

        angle_param_dict = {tuple(included_angle_type[0]): list(param[9:12])}

        predicted_f = []
        my_dict = copy.deepcopy(info_dict)
        for config_tag, info_list in my_dict.items():
            ff_forces = get_FF_forces(info_list[0], info_list[1], info_list[2], info_list[3], bond_param_dict,
                                      info_list[4], info_list[5], angle_param_dict)[10:13]
            predicted_f.append([force_list for force_list in ff_forces])
        print(predicted_f)
        residue = np.reshape(np.array(np.reshape(predicted_f, [-1, 3])) - np.array(np.reshape(DFT_f_on_O, [-1, 3])), -1)
        print(np.mean(residue ** 2))
        return np.mean(residue ** 2) * 100

    def get_fitting_parameters(initial_param, info_dict):
        # initial_param = np.array(list(initial_param_dict.values())).reshape(-1)
        res = minimize(get_residue, initial_param, method='Powell', args=info_dict, options={'ftol': 0.13})
        # res = least_squares(get_residue, initial_param)  # bounds=(0, np.inf)
        print(res.success)
        return res.x

    my_dict = copy.deepcopy(info_dict)
    print(get_fitting_parameters([*param1, *param2], my_dict))
    # does not converge, try including the angular term
    # print(get_residue([1.0, 1.0, 0.1, 1.1, 2.0, 0.1, 1.2, 4.0, 0.2, 2.3, 100.0]))
    """
