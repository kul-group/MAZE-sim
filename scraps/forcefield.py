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
from statistics import mode
import pickle
import time


def get_EF_atom_indices(atoms):
    """
    for index tracking, to ensure we are comparing the DFT and FF forces on the same EF atoms after before and after
    scooping out the smaller cluster.
    alse used for recentering the cluster based on the EF-O atom
    """
    TM_list = ['Pt', 'Cu', 'Co', 'Pd', 'Fe', 'Cr', 'Rh', 'Ru']
    index_EF_TM = [a.index for a in atoms if a.symbol in TM_list]
    index_Al = [a.index for a in atoms if a.symbol == 'Al']
    nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
    nl.update(atoms)
    Al_neigh_list = np.concatenate((nl.get_neighbors(index_Al[0])[0], nl.get_neighbors(index_Al[1])[0]))
    Al_neigh_list = [x for x in Al_neigh_list if atoms[x].symbol == 'O']

    TM_neigh_list = np.concatenate((nl.get_neighbors(index_EF_TM[0])[0], nl.get_neighbors(index_EF_TM[1])[0]))
    centering_o = [[x for x in TM_neigh_list if list(TM_neigh_list).count(x) > 1 and x not in Al_neigh_list][0]]
    return index_EF_TM + centering_o


def get_capped_cluster(atoms, folder_path, file_name, save_traj, EF_O_index):
    """ #TODO: check whether capping is necessary
    Inconsistent capping (remove all caps for now, does not need this cluster to be physical)
    Possible fix: change mult in neighbor list

    Extract smaller cluster containing the extra-framework atoms and cap all the O. Then the capped cluster is moved
    to the center of the cell to avoid boundary issue.
    Save cluster in both .traj file and .pdb format.
    :param atoms:
    :param folder_path:
    :param file_name:
    :param save_traj: if True, save clusters into .traj as well, for later comparison and trouble shooting
    :param EF_O_index: if not none, will use this value, else, will find the index using Extraframework code
    :return: 1. EF-cluster including 13 atoms, index of the EF atoms in original zeolite, index of the EF atoms in
    the current cluster (the later two output index lists share the ordering)
    """
    EFMaker = ExtraFrameworkAnalyzer(atoms)
    cluster = atoms[[index for index in EFMaker.get_extraframework_cluster(EF_O_index)]]

    cluster_EF_index = get_EF_atom_indices(cluster)
    centering_pos = cluster.get_positions()[cluster_EF_index[-1]]
    recentered_cluster = EFMaker.recentering_atoms(cluster, centering_pos)[0]
    # FIXME: recentering doesn't work well for very small unit cells. eg. SOD
    # cluster = Zeolite(cluster).cap_atoms()

    proteindatabank.write_proteindatabank(folder_path + '/%s.pdb' % file_name, recentered_cluster)
    if save_traj is True:
        write(folder_path + '/%s.traj' % file_name, recentered_cluster)

    return cluster, EFMaker.get_extraframework_cluster(EF_O_index), cluster_EF_index


def label_pdb(folder_path, file_name, del_unlabeled_pdb):
    """
    Relabeling the Atom name in proteindatabank file. (required step for openMM)
    The same atom type connecting to different neighboring types are treated differently due to differences in their
    chemical environments, and is therefore named separately.
    :param folder_path:
    :param file_name:
    :param del_unlabeled_pdb:
    """
    filein = open(folder_path + '/%s.pdb' % file_name, 'r')
    fileout = open(folder_path + '/%s_labeled.pdb' % file_name, 'w')

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
    if del_unlabeled_pdb is True:
        os.remove(folder_path + '/%s.pdb' % file_name)


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
    class_O_EF = [get_EF_atom_indices(cluster)[-1]]
    class_O_Cu = [atom.index for atom in cluster if atom.symbol == 'O' and atom.index not in class_O_EF and
                  all(val not in class_H for val in nl.get_neighbors(atom.index)[0])]
    class_O_H = [atom.index for atom in cluster if atom.symbol == 'O' and atom.index not in class_O_EF + class_O_Cu]

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


def set_up_openMM_system(folder_path, cluster_tag_number, shortened_bond_list):
    """ Feed pdb topology file and xml force field file into openMM, generate a system for the MD simulation/force
    calculation.
    :param folder_path:
    :param cluster_tag_number:
    :param shortened_bond_list:
    :return pdb:
    :return system:
    """
    pdb = PDBFile(folder_path + '/cluster_%s_labeled.pdb' % cluster_tag_number)
    atoms = list(pdb.topology.atoms())

    for index in shortened_bond_list:
        pdb.topology.addBond(atoms[index[0]], atoms[index[1]])
    bonds = list(pdb.topology.bonds())

    write_xml(atoms, bonds, folder_path + '/forcefield.xml')
    FF = ForceField(folder_path + '/forcefield.xml')
    system = FF.createSystem(pdb.topology)
    return pdb, system


def custom_openMM_force_object(system, bond_list, bond_type_index_dict, bond_param_dict, angle_list=None,
                               angle_type_index_dict=None, angle_param_dict=None):
    """ #todo: add argument allowing this custom function to be fed in as an input (more flexible used-designed ff)
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
                try:
                    force.addBond(int(bond[0]), int(bond[1]), bond_param_dict.get(my_type))
                except:
                    my_type = tuple(reversed(my_type))
                    force.addBond(int(bond[0]), int(bond[1]), bond_param_dict.get(my_type))
                    # note: consider updating the info_dict to make it order insensitive
    system.addForce(force)

    force = HarmonicAngleForce()  # Harmonic angle
    for angle in angle_list:
        for my_type, my_index in angle_type_index_dict.items():
            if any(list(val) in my_index for val in list(permutations(angle))):
                type_tag = [tuple(val) for val in list(angle_param_dict.keys()) if val in list(permutations(my_type))]
                force.addAngle(int(angle[0]), int(angle[1]), int(angle[2]), *angle_param_dict.get(type_tag[0]))
    system.addForce(force)

    return system


def get_openMM_forces(pdb, system, bond_list, bond_type_index_dict, bond_param_dict, angle_list=None,
                      angle_type_index_dict=None, angle_param_dict=None):
    """ forces for a single configuration
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


# NOTE: section below deals with multiple input structures for force field training

def get_EF_O_index(traj):
    """
    get the mode of EF_O, and use that to extract the EF cluster for the force field training
    all EF atoms should have the same indices regardless of there is binds on the zeolite, as long as the zeolite
    framework is the same - (all EF atoms, aka. Cu-O-Cu insertion follows the same procedures)
    :param traj: traj of configurations containing all atoms, including both the zeolite backbone and EF atoms
    """
    EF_O_index_list = []
    for atoms in traj:
        try:
            EFAnalyzer = ExtraFrameworkAnalyzer(atoms)
            EF_O_index_list.append(EFAnalyzer.get_extraframework_cluster()[-1])
        except:
            ...
    return mode(tuple(EF_O_index_list))


def prep_topologies(folder_path, sample_zeolite, traj_name=None, save_traj=False, del_unlabeled_pdb=False,
                    show_all=False):
    """
    :param folder_path:
    :param sample_zeolite:
    :param traj_name:
    :param save_traj:
    :param del_unlabeled_pdb:
    :param show_all:
    """
    if traj_name is not None:
        traj = read(folder_path + '/%s.traj' % traj_name, ':')
        output_dir = os.path.join(folder_path, traj_name)
    else:
        traj = read(folder_path + '/%s.traj' % sample_zeolite, ':')
        output_dir = os.path.join(folder_path, sample_zeolite)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    cluster_traj, EF_O_index, EF_atoms_index, cluster_EF_index = [], get_EF_O_index(traj[0:100]), [], []
    for count, atoms in enumerate(traj):
        cluster, EF_atoms_index, cluster_EF_index = get_capped_cluster(atoms, output_dir, 'cluster_' + str(count),
                                                                       save_traj, [EF_O_index])
        cluster_traj.append(cluster)
        print(sample_zeolite, count)

    if show_all is True:
        view(cluster_traj)
    for val in range(len(traj)):
        label_pdb(output_dir, 'cluster_%s' % str(val), del_unlabeled_pdb)

    return EF_atoms_index, cluster_EF_index


def reformat_inputs(bond_param_dict, angle_param_dict):
    """ reformat input dict into lists
    :return bond_type: List[List[str]] eg. ['Cu', 'O']
    :return angle_type: List[List[str]] eg. ['Cu', 'O', 'Cu']
    :return param_list: List[float], extend all parameters into a single list, since scipy.optimize.minimize can only
    take an 1D array as initial guess parameter
    """
    bond_type, angle_type, param_list = [], [], []
    for types, indices in bond_param_dict.items():
        bond_type.append(list(types))
        param_list.extend([val for val in np.array(indices)])

    for types, indices in angle_param_dict.items():
        angle_type.append(list(types))
        param_list.extend([val for val in np.array(indices)])

    return bond_type, angle_type, param_list


def get_required_objects_for_ff(folder_path, cluster_tag_number, included_bond_type, included_angle_type,
                                print_numb_unique_bonds=False, print_numb_unique_angles=False):
    """ To reduce computational cost, objects such as pdb, system, shortened_bond_list, bond_type_index_dict are kept
    fixed for each configuration during the optimization (only run once).
    """
    cluster = read(folder_path + '/cluster_%s_labeled.pdb' % cluster_tag_number, '0')

    bond_index_list, shortened_bond_index_list = get_bonds(cluster, mult=2)
    bond_type_dict, whole_bond_type_list = get_property_types(cluster, bond_index_list)
    if print_numb_unique_bonds is True:
        print('Number of unique bond types:', len(bond_type_dict))

    angle_index_list, shortened_angle_index_list = get_angles(cluster, mult=2)
    angle_type_dict, whole_angle_type_list = get_property_types(cluster, angle_index_list)
    if print_numb_unique_angles is True:
        print('Number of unique angle types:', len(angle_type_dict))

    bond_type_index_dict = get_type_index_pair(bond_type_dict, whole_bond_type_list, bond_index_list)
    angle_type_index_dict = get_type_index_pair(angle_type_dict, whole_angle_type_list, angle_index_list)
    shortened_bond_list = shorten_index_list_by_types(bond_type_index_dict, include_property_type=included_bond_type)
    shortened_angle_list = shorten_index_list_by_types(angle_type_index_dict, include_property_type=included_angle_type)
    pdb, system = set_up_openMM_system(folder_path, cluster_tag_number, shortened_bond_list)

    return pdb, system, shortened_bond_list, bond_type_index_dict, shortened_angle_list, angle_type_index_dict


def get_FF_forces(param, info_dict, ini_bond_param_dict, ini_angle_param_dict, EF_index):
    """ openMM forces for multiple configuration based on the same set of parameters
    """
    bond_param_dict, angle_param_dict, number_of_bond_param = {}, {}, 0
    for count, (types, indices) in enumerate(ini_bond_param_dict.items()):
        bond_param_dict[types] = list(param[count * len(indices):(count + 1) * len(indices)])
        number_of_bond_param += len(indices)

    for count, (types, indices) in enumerate(ini_angle_param_dict.items()):
        angle_param_dict[types] = list(
            param[count * len(indices) + number_of_bond_param:(count + 1) * len(indices) + number_of_bond_param])

    predicted_f = []
    my_dict = copy.deepcopy(info_dict)
    for config_tag, info_list in my_dict.items():
        ff_forces = get_openMM_forces(info_list[0], info_list[1], info_list[2], info_list[3], bond_param_dict,
                                      info_list[4], info_list[5], angle_param_dict)[EF_index]  # fixme: EF_index hardcoded
        predicted_f.append([force_list for force_list in ff_forces])

    return predicted_f


def get_DFT_forces_single(atoms, atom_index):
    """
    reference DFT forces on single atoms
    """
    f_vec = atoms.calc.results['forces'][atom_index]  # self.atoms.get_forces()[atom_index]
    f_mag = np.linalg.norm(f_vec)
    return f_vec


def get_residue(param, info_dict, DFT_f, ini_bond_param_dict, ini_angle_param_dict, EF_index):
    """optimize force field parameters by minimizing this function
    """
    predicted_f = get_FF_forces(param, info_dict, ini_bond_param_dict, ini_angle_param_dict, EF_index)
    residue = np.reshape(np.array(np.reshape(predicted_f, [-1, 3])) - np.array(np.reshape(DFT_f, [-1, 3])), -1)
    print(np.mean(residue ** 2))
    return np.mean(residue ** 2)


def get_fitting_parameters(initial_param, info_dict, DFT_f, ini_bond_param_dict, ini_angle_param_dict, EF_index):
    """ # todo: consider adding constrains on the parameters, eg. Harmonic angle below np.pi,
    # todo: need a smart way to save all bounds
    """
    res = minimize(get_residue, initial_param, method='Powell', options={'ftol': 0.01, 'maxiter': 1000},
                   args=(info_dict, DFT_f, ini_bond_param_dict, ini_angle_param_dict, EF_index))
    print(res.success)
    return res


def make_parity_plot(ff_forces, dft_forces, atom_name):
    """ plot FF forces vs. DFT forces
    """
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
    plt.title('Force fitting on %s' % atom_name, fontsize=18)
    plt.show()


def temp_func():
    """
    # topologies prep
    EF_index_dict, cluster_EF_index_dict = {}, {}
    zeo_list = ['CHA', 'AEI', 'RHO', 'MWW', 'BEA', 'LTA', 'MAZ', 'MFI', 'MOR', 'SOD']
    for zeolite in zeo_list[0:3]:
        folder_path, sample_zeolite, traj_name = '/Users/jiaweiguo/Box/openMM_FF', zeolite, zeolite + '_minE'
        EF_index, cluster_EF_index = prep_topologies(folder_path, sample_zeolite, traj_name, del_unlabeled_pdb=True)
        EF_index_dict[zeolite] = EF_index
        cluster_EF_index_dict[zeolite] = cluster_EF_index

    # write all original EF-atom indices into dict for later extraction of the DFT forces
    with open('/Users/jiaweiguo/Box/openMM_FF/EF_index_dict.pickle', 'wb') as f:
        pickle.dump(EF_index_dict, f)

    # write all EF-atom indeices in smaller cluster into dict for later extraction of the FF forces
    with open('/Users/jiaweiguo/Box/openMM_FF/cluster_EF_index_dict.pickle', 'wb') as f:
        pickle.dump(cluster_EF_index_dict, f)
    """

    # input parameters prep
    ini_bond_param_dict = {('O-Cu', 'Cu'): [1.2, 4, 0.2], ('O-EF', 'Cu'): [1.2, 4, 0.2]}  # ('Al', 'Cu'): [1.2, 4, 0.2],
    ini_angle_param_dict = {('Cu', 'O-EF', 'Cu'): [2.3, 40], ('Cu', 'O-Cu', 'Cu'): [2.3, 40]}
    # angle:radians, k:kJ/mol/radian^2

    included_bond_type, included_angle_type, ini_param = reformat_inputs(ini_bond_param_dict, ini_angle_param_dict)
    print(included_bond_type, included_angle_type, ini_param)

    # save openMM properties of each configuration into dict
    zeolite = 'CHA'
    folder_path, sample_zeolite, traj_name = '/Users/jiaweiguo/Box/openMM_FF', zeolite, zeolite + '_minE'
    traj = read(folder_path + '/%s.traj' % traj_name, ':')

    """
    info_dict, output_path = {}, os.path.join(folder_path, traj_name)
    for cluster_tag_number in range(len(traj)):
        pdb, system, shortened_bond_list, bond_type_index_dict, shortened_angle_list, angle_type_index_dict = \
            get_required_objects_for_ff(output_path, cluster_tag_number, included_bond_type, included_angle_type)
        info_dict[cluster_tag_number] = [pdb, system, shortened_bond_list, bond_type_index_dict, shortened_angle_list,
                                         angle_type_index_dict]

    with open(output_path+'/info_dict.pickle', 'wb') as f:
        pickle.dump(info_dict, f)
    """

    with open(folder_path + '/EF_index_dict.pickle', 'rb') as f:
        EF_index_dict = pickle.load(f)

    DFT_f = []
    for atoms in traj:
        DFT_f.append([get_DFT_forces_single(atoms, atom_index=val) for val in EF_index_dict.get(zeolite)[-3:]])
    print(np.array(DFT_f).shape)

    with open(os.path.join(folder_path, traj_name) + '/info_dict.pickle', 'rb') as f:
        info_dict = pickle.load(f)

    with open(folder_path + '/cluster_EF_index_dict.pickle', 'rb') as f:
        cluster_EF_index_dict = pickle.load(f)

    # my_dict = copy.deepcopy(info_dict)
    # print(get_FF_forces(ini_param, my_dict, ini_bond_param_dict, ini_angle_param_dict, cluster_EF_index_dict.get(zeolite)))

    my_dict = copy.deepcopy(info_dict)  # important, need to keep openMM "systems" fixed
    res = get_fitting_parameters(ini_param, my_dict, DFT_f, ini_bond_param_dict, ini_angle_param_dict,
                                 cluster_EF_index_dict.get(zeolite))

    print([np.around(float(val), decimals=3) for val in res.x])
    FF_f = get_FF_forces(res.x, info_dict, ini_bond_param_dict, ini_angle_param_dict,
                         cluster_EF_index_dict.get(zeolite))
    make_parity_plot(np.array(np.reshape(FF_f, [-1, 3])), np.array(np.reshape(DFT_f, [-1, 3])), 'Cu-O-Cu')

    # pretty_print(angle_type_index_dict)


def initial_traj_trimming():
    # topologies prep
    EF_index_dict, cluster_EF_index_dict = {}, {}
    zeo_list = ['CHA', 'AEI', 'RHO', 'MWW', 'BEA', 'LTA', 'MAZ', 'MFI', 'MOR', 'SOD']
    for zeolite in zeo_list:
        folder_path, sample_zeolite, traj_name = '/Users/jiaweiguo/Box/openMM_FF', zeolite, zeolite + '_minE'
        EF_index, cluster_EF_index = prep_topologies(folder_path, sample_zeolite, traj_name, del_unlabeled_pdb=True)
        EF_index_dict[zeolite] = EF_index
        cluster_EF_index_dict[zeolite] = cluster_EF_index

    # write all original EF-atom indices into dict for later extraction of the DFT forces
    with open('/Users/jiaweiguo/Box/openMM_FF/EF_index_dict.pickle', 'wb') as f:
        pickle.dump(EF_index_dict, f)

    # write all EF-atom indeices in smaller cluster into dict for later extraction of the FF forces
    with open('/Users/jiaweiguo/Box/openMM_FF/cluster_EF_index_dict.pickle', 'wb') as f:
        pickle.dump(cluster_EF_index_dict, f)

    with open('/Users/jiaweiguo/Box/openMM_FF/EF_index_dict.pickle', 'rb') as f:
        EF_index_dict = pickle.load(f)

    # fixme: pdb relabeling scripts failed for 'MWW'

    all_traj = []
    zeo_list = ['CHA', 'AEI', 'RHO', 'MWW', 'BEA', 'LTA', 'MAZ', 'MFI', 'MOR', 'SOD']
    for zeolite in zeo_list:
        DFT_f = []
        traj = read('/Users/jiaweiguo/Box/openMM_FF/%s_all.traj' % zeolite, ":")
        for atoms in traj:
            DFT_f.append([get_DFT_forces_single(atoms, atom_index=val) for val in EF_index_dict.get(zeolite)[-3:]])
        round_num, f_lim, included_list, config_tag_list = 1, 4, [[], [], []], []
        for count, force_list in enumerate(DFT_f):
            for force_vec in force_list:
                if all([np.round(force_vec[index], round_num) not in included_list[index] for index in range(3)]) and \
                        all([abs(force_vec[index]) < f_lim for index in range(3)]):
                    [included_list[index].append(np.round(force_vec[index], round_num)) for index in range(3)]
                    config_tag_list.append(count)

        print(zeolite, len(list(np.unique(config_tag_list))))
        [all_traj.append(traj[index]) for index in list(np.unique(config_tag_list))]

    write('/Users/jiaweiguo/Box/openMM_FF/reduced_all.traj', all_traj)

    # plt.scatter(DFT_f, DFT_f)
    # plt.show()


if __name__ == '__main__':
    # construct the bond_type_index_dict and angle_type_index_dict beforehand, all cluster will share one, expand if 
    # run into new types
    
    
    
    
    """
    tic = time.perf_counter()
    zeolite = 'CHA'
    folder_path, sample_zeolite, traj_name = '/Users/jiaweiguo/Box/openMM_FF', zeolite, zeolite + '_md'
    ini_bond_param_dict = {('O-Cu', 'Cu'): [1.2, 4, 0.2], ('O-EF', 'Cu'): [1.2, 4, 0.2]}  # ('Al', 'Cu'): [1.2, 4, 0.2],
    ini_angle_param_dict = {('Cu', 'O-EF', 'Cu'): [2.3, 40], ('Cu', 'O-Cu', 'Cu'): [2.3, 40]}
    included_bond_type, included_angle_type, ini_param = reformat_inputs(ini_bond_param_dict, ini_angle_param_dict)

    info_dict, output_path = {}, os.path.join(folder_path, traj_name)
    files = [files for files in os.listdir(os.path.join(folder_path, traj_name)) if '.pdb' in files]
    for cluster_tag_number in range(5):  # len(files)
        pdb, system, shortened_bond_list, bond_type_index_dict, shortened_angle_list, angle_type_index_dict = \
            get_required_objects_for_ff(output_path, cluster_tag_number, included_bond_type, included_angle_type)
        info_dict[cluster_tag_number] = [pdb, system, shortened_bond_list, bond_type_index_dict, shortened_angle_list,
                                         angle_type_index_dict]
        # print(cluster_tag_number)
    toc = time.perf_counter()
    print(f"Program terminated in {toc - tic:0.4f} seconds")
    """
    
