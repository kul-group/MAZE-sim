from maze.zeotypes import Zeotype, ImperfectZeotype
from ase import Atoms
from typing import Union, Tuple
from collections import defaultdict
from ase.neighborlist import natural_cutoffs, NeighborList
from ase import Atoms
import numpy as np
from ase.visualize import view
import copy as copy


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
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None, cif_dir=None):
        self.zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
        self.EFzeolite = ImperfectZeotype(self.zeolite)
        self.traj_1Al = []
        self.traj_2Al = []
        self.count = 0

    def get_unique_t_sites(self) -> dict:
        # return a dictionary containing all unique T site names and indices
        unique_t_site_indices = {}
        for site_name, value in self.zeolite._site_to_atom_indices.items():
            if 'T' in site_name:
                unique_t_site_indices[site_name] = value[1]
        return unique_t_site_indices

    def create_1Al_replacement(self):
        # create structures with 1 Al replacement for all unique T sites
        unique_t_site_indices = self.get_unique_t_sites()
        for site_name, t_site in unique_t_site_indices.items():
            position = self.EFzeolite.get_positions()[t_site]
            new_zeolite = self.EFzeolite.delete_atoms([t_site])
            new_zeolite = new_zeolite.add_atoms(Atoms('Al', positions=[position]), 'Al')
            self.traj_1Al.append(new_zeolite)

    def create_2Al_replacement(self):
        # make the second Al replacement
        done_indices = []
        for atoms in self.traj_1Al:
            ini_atoms = atoms
            index_Al = [a.index for a in atoms if a.symbol == 'Al']

            # skip process if there are already 2 Al replacement
            if len(index_Al) != 1:
                print('Error: No more than 1 Al in the structure.')
                atoms[index_Al[len(index_Al) - 1]].symbol = 'Si'
                index_Al = [a.index for a in atoms if a.symbol == 'Al']

            T_index = [a.index for a in atoms if a.symbol in ['Al', 'Si']]
            for index in np.arange(len(T_index)):
                atoms = copy.deepcopy(ini_atoms)
                indices, offsets = atoms.neighbor_list.get_neighbors(T_index[index])
                if len(indices) != 4:
                    continue
                else:
                    if T_index[index] != index_Al and [T_index[index], index_Al] not in done_indices:
                        if 3.3 < atoms.get_distance(index_Al, T_index[index]) < 9.0:
                            atoms[T_index[index]].symbol = 'Al'
                            self.traj_2Al.append(atoms)
                            self.count += 1
                            done_indices.append([index_Al, T_index[index]])
                            done_indices.append([T_index[index], index_Al])

    @staticmethod
    def get_close_framework_atom_positions(centering_index, atoms, radius_cutoff):
        close_framework_atom_positions = []
        for count in range(len(atoms)):
            dist = atoms.get_distance(centering_index, count, mic=True)
            if dist <= radius_cutoff:
                close_framework_atom_positions.append(atoms[count].position)
        return close_framework_atom_positions

    def _get_reference_with_S_in_between_Al_pairs(self, atoms):
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

