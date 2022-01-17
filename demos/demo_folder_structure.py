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
            write(my_path + '.vasp', atoms_2Al, format='vasp')
            #write(my_path + '.traj', atoms_2Al)

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

    def func(atoms, EF_atoms, output_dir0, file, ref_list=None, ref_index=None, skip_rotation=False, min_cutoff=0,
             max_cutoff=6, zeolite_dist_cutoff=1.5, AlAl_dist_cutoff=9):
        EFzeolite = ExtraFrameworkMaker()
        my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(atoms, EF_atoms, ref_list, ref_index, skip_rotation, min_cutoff,
                                                        max_cutoff, zeolite_dist_cutoff, AlAl_dist_cutoff)
        count = 1
        if 1 < len(my_atoms) < 10:
            for each_atoms in my_atoms:
                output_dir1 = os.path.join(output_dir0, file + '_%s' % str(count))
                Path(output_dir1).mkdir(parents=True, exist_ok=True)
                write(output_dir1 + '/CuOCu.traj', each_atoms)
                count += 1
        else:
            output_dir1 = os.path.join(output_dir0, file)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            write(output_dir1 + '/CuOCu.traj', my_atoms)

    finished_file, failed_file = [], []
    sample_zeolite = 'CHA'
    filepath = '/Users/gjw123/Desktop/2Al_mic_True/00_%s/' % sample_zeolite
    output_dir0 = '/Users/gjw123/Desktop/2Al_CuOCu_mic_True/rerun_%s/' % sample_zeolite
    # trouble_list = ['T4_T12_216_280']
    files = [files for files in os.listdir(filepath) if 'T' in files and 'pickle' not in files]
    # in trouble_list and 'pickle' not in files]
    Path(output_dir0).mkdir(parents=True, exist_ok=True)
    for file in files:
        atoms = read(os.path.join(filepath, file, '%s.traj' % file), '0')
        # view(atoms)
        EF_atoms = read('/Users/gjw123/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
        EF_atoms.set_cell(atoms.get_cell())
        try:
            func(copy.deepcopy(atoms), copy.deepcopy(EF_atoms), output_dir0, file, max_cutoff=4)
        except:
            try:
                func(copy.deepcopy(atoms), copy.deepcopy(EF_atoms), output_dir0, file, max_cutoff=4, zeolite_dist_cutoff=1.2)
            except:
                try:
                    EF_atoms = read('/Users/gjw123/Desktop/MAZE-sim-master/demos/CuOCu_cluster_smaller.traj', '0')
                    EF_atoms.set_cell(atoms.get_cell())
                    func(copy.deepcopy(atoms), copy.deepcopy(EF_atoms), output_dir0, file, max_cutoff=5)
                except:
                    try:
                        func(copy.deepcopy(atoms), copy.deepcopy(EF_atoms), output_dir0, file, max_cutoff=5, zeolite_dist_cutoff=1.2)
                    except:
                        print(file, 'failed!')
                        failed_file.append(file)
        finished_file.append(file)
        print(file, 'is done!')

    """
    traj, finished_file, failed_file = [], [], []
    sample_zeolite = 'MFI'
    filepath = '/Users/jiaweiguo/Box/all_zeo_database/2Al_mic_True/00_%s/' % sample_zeolite
    output_dir0 = '/Users/jiaweiguo/Box/P1_pair_site/MFI_2Al_Pt_mic_True/'
    original_dir = '/Users/jiaweiguo/Box/P1_pair_site/MFI_2Al_Pt'
    original_files = [files for files in os.listdir(original_dir) if 'T' in files]
    files = [files for files in os.listdir(filepath) if 'T' in files]
    # files = [files for files in os.listdir(filepath) if files in ['T5_T9_224_258', 'T8_T11_248_274']]
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
        EF_atoms = Atoms('Pt', positions=[[5, 5, 5]])
        EF_atoms.set_cell(atoms.get_cell())
        # need to reconsider Al-Al distance cutoff, within a huge ring, won't be stabilized

        try:
            EFzeolite = ExtraFrameworkMaker()
            new_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms), skip_rotation=True,
                                                             zeolite_dist_cutoff=1, max_cutoff=1)
            output_dir1 = os.path.join(output_dir0, file)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            write(output_dir1 + '/starting.traj', new_atoms)
            finished_file.append(file)
            traj.append(new_atoms)
        except:
            try:
                EFzeolite = ExtraFrameworkMaker()
                new_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms), skip_rotation=True,
                                                                 zeolite_dist_cutoff=1, max_cutoff=4)
                output_dir1 = os.path.join(output_dir0, file)
                Path(output_dir1).mkdir(parents=True, exist_ok=True)
                write(output_dir1 + '/starting.traj', new_atoms)
                finished_file.append(file)
                traj.append(new_atoms)
            except:
                try:
                    EFzeolite = ExtraFrameworkMaker()
                    new_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms), skip_rotation=True,
                                                                    zeolite_dist_cutoff=1, max_cutoff=3)
                    output_dir1 = os.path.join(output_dir0, file)
                    Path(output_dir1).mkdir(parents=True, exist_ok=True)
                    write(output_dir1 + '/starting.traj', new_atoms)
                    traj.append(new_atoms)
                    finished_file.append(file)
                except:
                    failed_file.append(file)
                    print(file, 'failed!')
    """

    """
    traj, finished_file, failed_file = [], [], []
    sample_zeolite = 'MOR'
    filepath = '/Users/gjw123/Desktop/2Al_mic_True/00_%s/' % sample_zeolite
    output_dir0 = '/Users/gjw123/Desktop/2Al_Cu_zeo/00_%s/' % sample_zeolite
    #files = [files for files in os.listdir(filepath) if 'T' in files]
    files = [files for files in os.listdir(filepath) if files in ['T2_T3_112_128', 'T2_T2_112_119', 'T1_T3_96_128', 'T1_T1_96_103', 'T2_T3_112_132', 'T2_T4_112_136', 'T1_T1_96_105', 'T4_T4_136_143', 'T1_T2_96_116', 'T3_T3_128_135']]
    Path(output_dir0).mkdir(parents=True, exist_ok=True)
    for file in files:
        output_dir1 = os.path.join(output_dir0, file)
        Path(output_dir1).mkdir(parents=True, exist_ok=True)
        atoms = read(os.path.join(filepath, file, '%s.traj' % file), '0')
        # view(atoms)
        EF_atoms = Atoms('Cu', positions=[[5, 5, 5]])
        EF_atoms.set_cell(atoms.get_cell())
        try:
            EFzeolite = ExtraFrameworkMaker()
            new_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms), skip_rotation=True,
                                                             zeolite_dist_cutoff=1.5, max_cutoff=4.5)
            output_dir1 = os.path.join(output_dir0, file)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            write(output_dir1 + '/starting.traj', new_atoms)
            finished_file.append(file)
            traj.append(new_atoms)
        except:
            try:
                EFzeolite = ExtraFrameworkMaker()
                new_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms), skip_rotation=True,
                                                                 zeolite_dist_cutoff=1.5, max_cutoff=2)
                output_dir1 = os.path.join(output_dir0, file)
                Path(output_dir1).mkdir(parents=True, exist_ok=True)
                write(output_dir1 + '/starting.traj', new_atoms)
                finished_file.append(file)
                traj.append(new_atoms)
            except:
                try:
                    EFzeolite = ExtraFrameworkMaker()
                    new_atoms = EFzeolite.insert_ExtraFrameworkAtoms(copy.copy(atoms), copy.copy(EF_atoms), skip_rotation=True,
                                                                    zeolite_dist_cutoff=1.5, max_cutoff=3)
                    output_dir1 = os.path.join(output_dir0, file)
                    Path(output_dir1).mkdir(parents=True, exist_ok=True)
                    write(output_dir1 + '/starting.traj', new_atoms)
                    traj.append(new_atoms)
                    finished_file.append(file)
                except:
                    failed_file.append(file)
                    print(file, 'failed!')
    """

    """
    sample_zeolite = 'CHA'
    zeo_dir = '/Users/jiaweiguo/Box/all_zeo_database/Silica_unopt/00_%s/%s.traj' % (sample_zeolite, sample_zeolite)
    output_1Al = '/Users/jiaweiguo/Desktop/SSZ_13/1Al'
    output_2Al = '/Users/jiaweiguo/Desktop/SSZ_13/2Al/without_edges'
    # output_1Al = '/Users/jiaweiguo/Box/all_zeo_database/1Al/00_%s' % sample_zeolite
    # output_2Al = '/Users/jiaweiguo/Box/all_zeo_database/2Al_mic_True/00_%s' % sample_zeolite
    traj_saver = SaveExtraFrameworkConfigurations(sample_zeolite)
    traj_saver.save_all_Al_zeo(zeo_dir, output_1Al, output_2Al)
    """
