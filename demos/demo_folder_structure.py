# generating folder structures for all zeolites

from maze.extra_framework_maker import ExtraFrameworkMaker
from maze.io_zeolite import read_vasp
from maze.zeolite import PerfectZeolite, Zeolite
from ase.neighborlist import natural_cutoffs, NeighborList
import os
from pathlib import Path
from ase.io import write, read
from ase.visualize import view
import copy
import shutil
from glob import glob
from ase.constraints import FixAtoms
from ase import Atoms
import numpy as np
import pickle


class SaveExtraFrameworkConfigurations(object):

    def __init__(self, sample_zeolite: str):
        self.sample_zeolite = sample_zeolite

    def save_bare_zeo(self, output_dir):
        EFzeolite = ExtraFrameworkMaker(iza_code=self.sample_zeolite)
        output_dir0 = os.path.join(output_dir, '00_' + sample_zeolite)
        Path(output_dir0).mkdir(parents=True, exist_ok=True)
        my_path = os.path.join(output_dir0, sample_zeolite)
        write(my_path + '.traj', EFzeolite.EFzeolite)
        print(self.sample_zeolite, 'is done!')

    def save_all_Al_zeo(self, zeo_dir, output_1Al, output_2Al):
        EFzeolite = ExtraFrameworkMaker(self.sample_zeolite, zeo_dir)
        EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
        info_dict = {'T_indices': EFzeolite.t_site_indices, 'Atoms_count': len(EFzeolite.EFzeolite),
                     'Al_pair_count': len(EFzeolite.traj_2Al), 'T_count': len(EFzeolite.t_site_indices)}
        print('Number of unique Al pairs for %s: ' % self.sample_zeolite, len(EFzeolite.traj_2Al))
        """
        print(len(EFzeolite.traj_1Al))
        if not os.path.exists(output_1Al):
            os.mkdir(output_1Al)
        for site_name, atoms_1Al in EFzeolite.dict_1Al_replaced.items():
            output_dir1 = os.path.join(output_1Al, site_name)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            my_path = os.path.join(output_dir1, site_name)
            write(my_path + '.traj', atoms_1Al[0])
        """
        if not os.path.exists(output_2Al):
            os.mkdir(output_2Al)
        for site_name, atoms_2Al in EFzeolite.dict_2Al_replaced.items():
            output_dir2 = os.path.join(output_2Al, site_name)
            Path(output_dir2).mkdir(parents=True, exist_ok=True)
            my_path = os.path.join(output_dir2, site_name)
            write(my_path + '.traj', atoms_2Al)

        with open(output_2Al + '/info_dict_%s.pickle' % self.sample_zeolite, 'wb') as f:
            pickle.dump(info_dict, f)

    def save_H_structures(self, dir_name, Al_num: str):
        folder = dir_name + '/' + self.sample_zeolite + '/' + Al_num
        filepaths = [os.path.join(folder, dirs) for dirs in os.listdir(folder_name) if 'T' in dirs]
        for filepath in filepaths:
            files = [files for files in os.listdir(filepath) if '.traj' in files]
            for file in files:
                atoms = read(os.path.join(filepath, file))
                EFzeolite = ExtraFrameworkMaker()
                dict_H = EFzeolite.get_Bronsted_sites(atoms)
                for site_name, my_atoms in dict_H.items():
                    output_dir1 = os.path.join(filepath, 'H')
                    Path(output_dir1).mkdir(parents=True, exist_ok=True)
                    my_path = os.path.join(output_dir1, self.sample_zeolite + '_' + Al_num + '_' + site_name)
                    write(my_path + '.traj', my_atoms)

    def save_all_ZCu(self, filepath, output_dir):
        filepath = filepath + '/00_' + self.sample_zeolite
        files = [files for files in os.listdir(filepath) if 'T' in files]
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        for file in files:
            atoms = read(os.path.join(filepath, file, '%s.traj' % file), '0')
            EFzeolite = ExtraFrameworkMaker()
            dict_ZCu = EFzeolite.get_Z_TM(atoms, 2.6, 'Cu')

            for site_name, my_atoms in dict_ZCu.items():
                output_dir1 = os.path.join(output_dir, file + '_' + site_name)
                Path(output_dir1).mkdir(parents=True, exist_ok=True)
                write(output_dir1 + '/1Cu.traj', my_atoms)
        print(self.sample_zeolite, 'is done!')

    def save_all_CuOCu(self, filepath, output_dir0):
        files = [files for files in os.listdir(filepath) if 'T' in files]
        Path(output_dir0).mkdir(parents=True, exist_ok=True)
        for file in files:
            output_dir1 = os.path.join(output_dir0, file)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            atoms = read(os.path.join(filepath, file, '%s.traj' % file), '0')
            EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
            EF_atoms.set_cell(atoms.get_cell())
            try:
                EFzeolite = ExtraFrameworkMaker()
                my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(atoms, EF_atoms, max_cutoff=6)
                write(output_dir1 + '/CuOCu.traj', my_atoms)
            except:
                print(self.sample_zeolite, file, ' failed!')


if __name__ == '__main__':
    topo = ['CHA', 'MOR', 'MWW', 'FER', 'MFI', 'BEC', 'MON', 'MSE', 'AFO', 'AHT', 'BOG', 'CFI', 'JOZ', 'JRY', 'JSN',
            'CGF', 'DON', 'EMT', 'EON', 'EUO', 'GON', 'IWS', 'LTF', 'MTW', 'OBW', 'OSI', 'RHO', 'RSN', 'JST', 'JSW',
            'SBE', 'SSY', 'TER', 'VFI', 'WEI', 'ABW', 'ACO', 'AET', 'AFI', 'AFN', 'AFR', 'AFT', 'AFY', 'ITW', 'IRR',
            'ANA', 'APC', 'APD', 'AST', 'ASV', 'ATN', 'ATO', 'ATT', 'ATV', 'AWW', 'BCT', 'BIK', 'BOF', 'BEA', 'LTJ',
            'BRE', 'BSV', 'CAN', 'CGS', 'CHI', 'CON', 'CZP', 'DDR', 'DFO', 'DFT', 'EAB', 'EDI', 'SAF', 'RTE', 'PUN',
            'ERI', 'ESV', 'ETR', 'EZT', 'FAR', 'FRA', 'GIS', 'GIU', 'GME', 'GOO', 'HEU', 'IFR', 'IHW', 'PCR', 'OWE',
            'ISV', 'ITE', 'ITH', 'ITR', 'IWR', 'JBW', 'KFI', 'LEV', 'LIO', 'LOS', 'LOV', 'LTA', 'PAR', 'NPT', 'MVY',
            'LTL', 'MAR', 'MAZ', 'MEL', 'MEP', 'MER', 'MFS', 'MOZ', 'MTF', 'MTN', 'NAB', 'NAT', 'MSO', 'MEI', 'LIT',
            'NES', 'NON', 'NPO', 'NSI', 'OFF', 'OSO', 'PHI', 'PON', 'RRO', 'RTH', 'RUT', 'RWR', 'LAU', 'UEI', 'UFI',
            'RWY', 'SAO', 'SAS', 'SAV', 'SBN', 'SBS', 'SFE', 'SFF', 'SFG', 'SFH', 'SFN', 'SFO', 'UOS', 'UOZ', 'SZR',
            'SFS', 'SIV', 'SOD', 'STF', 'STW', 'THO', 'TOL', 'USI', 'UTL', 'VET', 'VNI', 'VSV', 'WEN', 'STI', 'SVV',
            'YUG', 'ZON', 'IWW', 'STT', 'SVR', 'STO', 'AEI', 'AEL', 'AEN', 'AFG', 'AFS', 'AFX', 'SGT', 'SOF', 'SOS',
            'ATS', 'AVL', 'AWO', 'BOZ', 'BPH', 'CAS', 'CDO', 'DAC', 'DOH', 'EPI', 'UWY', 'TON', 'SSF', 'SAT']

    # ['FAU', 'CLO', 'IWV', 'LTN', 'PAU', 'SBT', 'TUN', 'TSC', 'ITV', 'IMF']
    """
    failed_zeo = []
    for sample_zeolite in ['STO', 'WEN']:
        try:
            traj_saver = SaveExtraFrameworkConfigurations(sample_zeolite)
            # traj_saver.save_bare_zeo('/Users/jiaweiguo/Box/00_bare_zeo/')
            input_dir = '/Users/jiaweiguo/Box/all_zeo_database/1Al'
            output_dir = '/Users/jiaweiguo/Box/all_zeo_database/1Al_Cu/00_%s' % sample_zeolite
            traj_saver.save_all_ZCu(input_dir, output_dir)
        except:
            failed_zeo.append(sample_zeolite)
    print(failed_zeo)
    """
    """
    failed_zeo = []
    for sample_zeolite in topo:
        try:
            zeo_dir = '/Users/jiaweiguo/Box/all_zeo_database/Silica_unopt/00_%s/%s.traj' % (sample_zeolite, sample_zeolite)
            output_1Al = '/Users/jiaweiguo/Box/all_zeo_database/1Al/00_%s' % sample_zeolite
            output_2Al = '/Users/jiaweiguo/Box/all_zeo_database/2Al/00_%s' % sample_zeolite
            traj_saver = SaveExtraFrameworkConfigurations(sample_zeolite)
            traj_saver.save_all_Al_zeo(zeo_dir, output_1Al, output_2Al)
            print(sample_zeolite + ' is done!')
        except:
            failed_zeo.append(sample_zeolite)
    print(failed_zeo)
    
    for sample_zeolite in ['RRO', 'PUN', 'MTW']:  # 'SFN', 'SFE', 'PON', 'SSY', 'SFF', 'CGF', 'SOF', 'RRO', 'PUN', 'MTW'
        filepath = '/Users/jiaweiguo/Box/all_zeo_database/2Al/00_%s/' % sample_zeolite
        output_dir = '/Users/jiaweiguo/Box/all_zeo_database/2Al_CuOCu/rerun/01_%s/' % sample_zeolite
        traj_saver = SaveExtraFrameworkConfigurations(sample_zeolite)
        traj_saver.save_all_CuOCu(filepath, output_dir)
        print(sample_zeolite, ' is done!')
    """

    traj, finished_file, failed_file = [], [], []
    sample_zeolite = 'AEI'
    filepath = '/Users/jiaweiguo/Box/all_zeo_database/2Al_mic_True/00_%s/' % sample_zeolite
    output_dir0 = '/Users/jiaweiguo/Box/all_zeo_database/2Al_CuOCu_mic_True/01_%s/' % sample_zeolite
    files = [files for files in os.listdir(filepath) if 'T' in files]
    Path(output_dir0).mkdir(parents=True, exist_ok=True)
    for file in files:
        output_dir1 = os.path.join(output_dir0, file)
        Path(output_dir1).mkdir(parents=True, exist_ok=True)
        atoms = read(os.path.join(filepath, file, '%s.traj' % file), '0')
        # view(atoms)
        EF_atoms = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
        # EF_atoms = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/CuOCu_cluster_smaller.traj', '0')
        EF_atoms.set_cell(atoms.get_cell())
        # EFzeolite = ExtraFrameworkMaker()
        # my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(atoms, EF_atoms, max_cutoff=4, AlAl_dist_cutoff=10)
        try:
            EFzeolite = ExtraFrameworkMaker()
            my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms))
            write(output_dir1 + '/CuOCu.traj', my_atoms)
            traj.append(my_atoms)
            finished_file.append(file)
        except:
            try:
                EFzeolite = ExtraFrameworkMaker()
                my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms),
                                                                zeolite_dist_cutoff=1)
                write(output_dir1 + '/CuOCu.traj', my_atoms)
                traj.append(my_atoms)
                finished_file.append(file)
            except:
                try:
                    EF_atoms = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/CuOCu_cluster_smaller.traj', '0')
                    EF_atoms.set_cell(atoms.get_cell())
                    EFzeolite = ExtraFrameworkMaker()
                    my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms))
                    write(output_dir1 + '/CuOCu.traj', my_atoms)
                    traj.append(my_atoms)
                    finished_file.append(file)
                except:
                    try:
                        EF_atoms = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/CuOCu_cluster_smaller.traj', '0')
                        EF_atoms.set_cell(atoms.get_cell())
                        EFzeolite = ExtraFrameworkMaker()
                        my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms), zeolite_dist_cutoff=1)
                        write(output_dir1 + '/CuOCu.traj', my_atoms)
                        traj.append(my_atoms)
                        finished_file.append(file)
                    except:
                        print(file, 'failed!')
                        failed_file.append(file)


    """
    traj, finished_file, failed_file = [], [], []
    sample_zeolite = 'MFI'
    filepath = '/Users/jiaweiguo/Box/all_zeo_database/2Al/00_%s/' % sample_zeolite
    output_dir0 = '/Users/jiaweiguo/Box/P1_pair_site/MFI_2Al_Pt/'
    trouble_list = ['T1_T8_192_252']
    files = [files for files in os.listdir(filepath) if files in trouble_list]
    Path(output_dir0).mkdir(parents=True, exist_ok=True)
    for file in files:
        output_dir1 = os.path.join(output_dir0, file)
        Path(output_dir1).mkdir(parents=True, exist_ok=True)
        atoms = read(os.path.join(filepath, file, '%s.traj' % file), '0')
        view(atoms)
        nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
        nl.update(atoms)
        index_Al = [a.index for a in atoms if a.symbol == 'Al']
        Al_neigh_list = np.concatenate((nl.get_neighbors(index_Al[0])[0], nl.get_neighbors(index_Al[1])[0]))
        Al_neigh_list = [x for x in Al_neigh_list if atoms[x].symbol == 'O']
        c = FixAtoms(indices=[atom.index for atom in atoms if atom.symbol != 'Al' and atom.index not in Al_neigh_list])
        atoms.set_constraint(c)
        my_atoms = copy.copy(atoms)
        # need to reconsider Al-Al distance cutoff, within a huge ring, won't be stabilized

        output_dir1 = os.path.join(output_dir0, file)
        Path(output_dir1).mkdir(parents=True, exist_ok=True)
        try:
            EF_atoms = Atoms('Pt', positions=[[0, 0, 0]])
            EF_atoms.set_cell(atoms.get_cell())
            EFzeolite = ExtraFrameworkMaker()
            new_atoms = EFzeolite.insert_ExtraFrameworkAtoms(my_atoms, EF_atoms, skip_rotation=True, max_cutoff=4)
            write(output_dir1 + '/starting.traj', new_atoms)
            finished_file.append(file)
            traj.append(new_atoms)
        except:
            failed_file.append(file)
            print(file, 'failed!')
    """

    """
    # consider adding a condition controling the distance of EF atoms from the Al, can not exceed certain criteria
    # or need to remove the occupancy check criteria
    traj = []
    filepath = '/Users/jiaweiguo/Box/P1_pair_site/MFI_2Al_Pt/'
    files = [files for files in os.listdir(filepath) if 'T' in files]
    for each_file in files:
        atoms = read(os.path.join(filepath, each_file, 'starting.traj'), '0')
        traj.append(atoms)
    view(traj)
    """

    """
    sample_zeolite = 'AEI'
    zeo_dir = '/Users/jiaweiguo/Box/all_zeo_database/Silica_unopt/00_%s/%s.traj' % (sample_zeolite, sample_zeolite)
    output_1Al = '/Users/jiaweiguo/Box/all_zeo_database/1Al/00_%s' % sample_zeolite
    output_2Al = '/Users/jiaweiguo/Box/all_zeo_database/2Al_mic_True/00_%s' % sample_zeolite
    traj_saver = SaveExtraFrameworkConfigurations(sample_zeolite)
    traj_saver.save_all_Al_zeo(zeo_dir, output_1Al, output_2Al)
    """

    """
    zeo_dir = '/Users/jiaweiguo/Box/all_zeo_database/Silica_unopt/00_%s/%s.traj' % (sample_zeolite, sample_zeolite)
    EFzeolite = ExtraFrameworkMaker(sample_zeolite, zeo_dir)
    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    info_dict = {'T_indices': EFzeolite.t_site_indices, 'Atoms_count': len(EFzeolite.EFzeolite),
                 'Al_pair_count': len(EFzeolite.traj_2Al), 'T_count': len(EFzeolite.t_site_indices)}
    print('Number of unique Al pairs for %s: ' % sample_zeolite, len(EFzeolite.traj_2Al))
    write('/Users/jiaweiguo/Box/zeo_dimer_formation_analysis/benchmarking_mic/mic_True.traj', EFzeolite.traj_2Al)
    """
