from maze.zeotypes import Zeotype, ImperfectZeotype
from ase import Atoms
from typing import Union, Tuple
from collections import defaultdict
from ase.neighborlist import natural_cutoffs, NeighborList
from ase import Atoms
import numpy as np
from ase.visualize import view
import copy as copy
from ase.io import write, read


'''
#TODO: At some point make a common ABC for the ExtraFrameworkAtoms and Adsorbate
class ExtraFrameworkAtoms(Atoms):
    """
    This class represents the extraframework atoms that will be added to the
    """
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None,  host_zeotype=None, name='', description=''):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions, cell,
                         pbc, celldisp, constraint, calculator, info, velocities)

        assert '_' not in description, 'cannot add _ to description'
        if isinstance(symbols, self.__class__):
            if host_zeotype is None:
                host_zeotype = symbols.host_zeotype
            if description == '':
                description = symbols.description
            if name == '':
                name = symbols.name

        self.host_zeotype = host_zeotype
        self.description = description
        self.name = name


class ExtraFramework(ImperfectZeotype):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None, cif_dir=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, site_to_atom_indices,
                         atom_indices_to_site, additions)


    def integrate_efa(self, efa: Atoms) -> Tuple['ExtraFramework', ExtraFrameworkAtoms]:
        """
        Integrate extraframework atoms into the extraframework and return new
        extraframework and the extraframework atoms
        :param efa:  Extraframework atoms needed for
        :type efa: Atoms, ExtraFrameworkAtoms
        :return:  Extraframework with efa integrated and ExtraFrameworkAtoms object
        :rtype: Tuple['ExtraFramework', ExtraFrameworkAtoms]
        """
        efa_name = 'extraframework_atoms'
        efa = ExtraFrameworkAtoms(efa)
        new_self = self.add_atoms(efa, efa_name, short_description=efa.description)
        efa.name = new_self.additions[efa_name][-1]
        efa.host_zeotype = new_self
        return new_self, efa
'''


class ExtraFramework(object):
    def __init__(self, cif_dir=None):
        """ #TODO: Add descriptions here
        :param cif_dir: cif file directory of the zeolite structures (IZA database)
        """
        if cif_dir is not None:
            self.EFzeolite = ImperfectZeotype.make(cif_dir)
        self.dict_t_sites_1Al_replaced = {}
        self.traj_1Al = []
        self.traj_2Al = []
        self.count_all_Al_pairs = 0
        self.TM_list = ['Pt', 'Cu', 'Co', 'Pd', 'Fe', 'Cr', 'Rh', 'Ru']

    ## need checking
    def make_extra_frameworks(self):
        self.replace_1Al()
        self.replace_2Al_unique_pairs()
        # add in Cu here

    def get_t_sites(self) -> dict:
        """ This function gets all unique T sites and corresponding atom indices for each T sites
        :return: dictionary mapping each unique T site tags with all atom indices with the same tag
        """
        t_site_indices = {}
        for site_name, value in self.EFzeolite.site_to_atom_indices.items():
            if 'T' in site_name:
                t_site_indices[site_name] = value
        return t_site_indices

    def replace_1Al(self):
        """ This function takes in a perfect zeolite and replace one Si atom with an Al for every T site
        :return: a dictionary of trajectories for each T site tags
        """
        t_site_indices = self.get_t_sites()
        for site_name, t_site in t_site_indices.items():
            traj_t_sites = []
            for index in t_site:
                pos = self.EFzeolite.get_positions()[index]
                new_zeo = self.EFzeolite.delete_atoms([index])
                new_zeo = new_zeo.add_atoms(Atoms('Al', positions=[pos]), 'Al')
                self.traj_1Al.append(new_zeo)
                traj_t_sites.append(new_zeo)
            self.dict_t_sites_1Al_replaced[site_name] = traj_t_sites

    def replace_2Al_unique_pairs(self, cutoff_radius=9):
        """ This function makes the 2 Al replacement for all possible pairs (not limited to unique T-site pairs since
        even though the binding energies might be the same, the geometric properties, such as, Al-Al distance, are
        different). Replacement obeys the Lowenstein's rule, and Al pairs with distance greater than the "cutoff_radius"
        are ignored.
        :param cutoff_radius: replace the second Si site within some cutoff radius (9 Angstrom by default) around the
        first replacement site which is done using function "replace_1Al".
        :return:
        """
        done_indices = []

        for zeolite in self.traj_1Al:
            atoms = ImperfectZeotype(zeolite)  # already have neighbor_list
            ini_atoms = copy.deepcopy(atoms)
            index_Al = [a.index for a in atoms if a.symbol == 'Al']
            '''
            # skip process if there are already 2 Al replacement
            if len(index_Al) != 1:
                print('Error: No more than 1 Al in the structure.')
                atoms[index_Al[len(index_Al) - 1]].symbol = 'Si'
                index_Al = [a.index for a in atoms if a.symbol == 'Al']
            '''
            neighboring_Si = []
            neigh_o_indices, offsets = atoms.neighbor_list.get_neighbors(index_Al[0])
            for each_neigh_o_index in neigh_o_indices:
                neigh_si_indices, offsets = atoms.neighbor_list.get_neighbors(each_neigh_o_index)
                for each_neigh_si_index in neigh_si_indices:
                    if each_neigh_si_index != index_Al[0]:
                        neighboring_Si.append(each_neigh_si_index)

            if len(neighboring_Si) != 4:
                continue
            else:
                Si_index = [a.index for a in atoms if a.symbol in ['Si']]
                for index in Si_index:
                    if index not in neighboring_Si and [index, index_Al] not in done_indices:
                        # first condition is the Lowenstein's rule
                        if 3.3 < atoms.get_distance(index_Al, index) < cutoff_radius:
                            atoms = copy.deepcopy(ini_atoms)
                            atoms[index].symbol = 'Al'
                            self.traj_2Al.append(atoms)
                            self.count_all_Al_pairs += 1
                            done_indices.append([index_Al, index])
                            done_indices.append([index, index_Al])

    @staticmethod
    def get_close_framework_atom_positions(centering_index, atoms, radius_cutoff):
        """
        :param centering_index:
        :param atoms:
        :param radius_cutoff:
        :return:
        """
        close_framework_atom_positions = []
        for count in range(len(atoms)):
            dist = atoms.get_distance(centering_index, count, mic=True)
            if dist <= radius_cutoff:
                close_framework_atom_positions.append(atoms[count].position)
        return close_framework_atom_positions

    def _get_reference_with_S_in_between_Al_pairs(self, atoms):
        """

        :param atoms:
        :return:
        """
        Al_index = [a.index for a in atoms if a.symbol in ['Al']]
        vec_Al = atoms.get_distance(Al_index[0], Al_index[1], mic=True, vector=True)
        mid_Al = vec_Al / 2 + atoms.get_positions()[Al_index[0]]
        atoms = atoms + Atoms('S', positions=[mid_Al])

        index_S = [a.index for a in atoms if a.symbol == 'S'][0]
        d_thres = 10.0  # might need to scale it with respect to the zeolite size

        close_framework_atom_positions = self.get_close_framework_atom_positions(index_S, atoms, d_thres)
        # print(close_framework_atom_positions)
        del atoms[[atom.index for atom in atoms if atom.symbol == 'S']]  # remove reference point

        max_iter = 100
        num_iter = 0
        closest_distance = 1.5  # radius of Si atom = 1.5 Ang
        u_Al = vec_Al / np.linalg.norm(vec_Al)  # unit direction along the Al-Al pair

        final_pos = []
        d_min = 0.0
        while num_iter <= max_iter:
            rand_vec = np.array([2.0 * np.random.random() - 1 for i in np.arange(3)])
            u_rand = rand_vec / np.linalg.norm(rand_vec)  # random unit vector
            u_dir = np.cross(u_Al, u_rand) / np.linalg.norm(
                np.cross(u_Al, u_rand))  # some random vector perpendicular to the Al-Al vector
            u_dir = int(2 * np.round(np.random.random()) - 1) * u_dir
            # This takes a random number between 0 and 1, rounds it (0 or 1) and then shifts it to (-1 or 1)

            step_size = 0.5 * np.random.random_sample() * d_thres  # Pick step size randomly from 0 to half radius
            trial_pos = np.array(mid_Al + u_dir * step_size)

            distance_from_close_atoms = (close_framework_atom_positions - trial_pos)
            d_min = min(np.linalg.norm(distance_from_close_atoms, axis=1))

            num_iter += 1
            if closest_distance <= d_min <= d_thres:
                final_pos = trial_pos
                break

        cell_dimension = np.matmul(np.ones(3), atoms.get_cell())
        mic_constrained_pos = []
        for index, item in enumerate(final_pos):
            if item > cell_dimension[index]:
                item -= cell_dimension[index]
            elif item < 0:
                item += cell_dimension[index]
            mic_constrained_pos.append(item)

        atoms = atoms + Atoms('S', positions=[mic_constrained_pos])
        return atoms

    # no longer needed
    def insert_CuOCu(self, atoms):
        atoms = self._get_reference_with_S_in_between_Al_pairs(atoms)

        index_S = [a.index for a in atoms if a.symbol == 'S'][0]
        pos_centering_O = atoms.get_positions()[index_S]

        del atoms[[atom.index for atom in atoms if atom.symbol == 'S']]
        atoms = atoms + Atoms('O', positions=[pos_centering_O])
        d_CuO = 2
        Al_index = [a.index for a in atoms if a.symbol in ['Al']]
        for count in range(len(Al_index)):
            dir_Cu = atoms.get_distance(index_S, Al_index[count], mic=True, vector=True)
            dir_Cu = dir_Cu / np.linalg.norm(dir_Cu)
            pos_Cu = dir_Cu * d_CuO + pos_centering_O
            atoms = atoms + Atoms('Cu', positions=[pos_Cu])
        return atoms

    # no longer needed 
    def insert_TM(self, atoms, metal_type):
        # find metal sites and insert in between Al and reference S atoms (which will be converted to O later)
        atoms = self._get_reference_with_S_in_between_Al_pairs(atoms)

        index_S = [a.index for a in atoms if a.symbol == 'S'][0]
        pos = atoms.get_positions()[index_S]

        del atoms[[atom.index for atom in atoms if atom.symbol == 'S']]
        atoms = atoms + Atoms(metal_type, positions=[pos])
        return atoms

    def insert_ExtraFrameworkAtoms(self, my_zeolite, EF_atoms):
        """ This function takes in a zeolite backbone and an extra-framework cluster with the same cell dimensions as
        the zeolite. First, move the cluster center-of-mass to the reference position (indicated using an S atom). If
        there are more than one TMs in the cluster, the cluster is rotated so that the TM-TM vector is aligned with the
        Al-Al vector. Last, insert the cluster into "my_zeolite".
        :param my_zeolite: zeolite backbone with 2 Al atoms close to each other
        :param EF_atoms: the extra-framework cluster to be inserted in between the Al pair
        :return:
        """
        atoms = self._get_reference_with_S_in_between_Al_pairs(my_zeolite)
        cell = atoms.get_cell()

        EF_center = EF_atoms.get_center_of_mass()
        index_S = [a.index for a in atoms if a.symbol == 'S'][0]
        ref_center = atoms.get_positions()[index_S]
        translate_vector = ref_center - EF_center
        EF_atoms.translate(translate_vector)
        EF_atoms.wrap()

        TM_index = [atom.index for atom in EF_atoms if atom.symbol in self.TM_list]
        Al_index = [atom.index for atom in atoms if atom.symbol == 'Al']

        if len(TM_index) == 2:
            ini_TM_vec = EF_atoms.get_positions()[TM_index[0]] - EF_atoms.get_positions()[TM_index[1]]
            final_vec = atoms.get_positions()[Al_index[0]] - atoms.get_positions()[Al_index[1]]
            EF_atoms.rotate(ini_TM_vec, final_vec, center=EF_center)
            # ini_TM_vec is rotated into final_vec

        pos = EF_atoms.get_positions()
        atom_name = ''.join(EF_atoms.get_chemical_symbols())
        atoms = atoms + Atoms(atom_name, positions=pos, cell=cell)

        del atoms[[atom.index for atom in atoms if atom.symbol == 'S']]
        return atoms


def main():
    zeolite = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/MFI_2Al_replaced.traj', '0')
    EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')

    EFzeolite = ExtraFramework()
    inserted_atoms = EFzeolite.insert_ExtraFrameworkAtoms(zeolite, EF_atoms)
    view(inserted_atoms)

    # write('MFI_Co.traj', inserted_atoms)


if __name__ == '__main__':
    main()
