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

    def save_all_Al_zeo(self, zeo_dir, output_1Al, output_2Al):
        EFzeolite = ExtraFrameworkMaker(self.sample_zeolite, zeo_dir)
        EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
        info_dict = {'T_indices': EFzeolite.t_site_indices, 'Atoms_count': len(EFzeolite.EFzeolite),
                     'Al_pair_count': len(EFzeolite.traj_2Al), 'T_count': len(EFzeolite.t_site_indices)}
        print('Number of unique Al pairs for %s: ' % self.sample_zeolite, len(EFzeolite.traj_2Al))

        print(len(EFzeolite.traj_1Al))
        if not os.path.exists(output_1Al):
            os.mkdir(output_1Al)
        for site_name, atoms_1Al in EFzeolite.dict_1Al_replaced.items():
            output_dir1 = os.path.join(output_1Al, site_name)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            my_path = os.path.join(output_dir1, site_name)
            write(my_path + '.traj', atoms_1Al[0])

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
        files = [files for files in os.listdir(filepath) if 'T' in files]
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        for file in files:
            output_dir1 = os.path.join(output_dir0, file)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            atoms = read(os.path.join(filepath, file, '%s.traj' % file), '0')
            EFzeolite = ExtraFrameworkMaker()
            dict_ZCu = EFzeolite.get_Z_TM(atoms, 2.6, 'Cu')

            for site_name, my_atoms in dict_ZCu.items():
                output_dir1 = os.path.join(output_dir, sub_dir_name + '_' + site_name)
                Path(output_dir1).mkdir(parents=True, exist_ok=True)
                write(output_dir1 + '/1Cu.traj', my_atoms)
        print(self.sample_zeolite, ' is done!')

    def save_all_CuOCu(self, filepath, output_dir0):
        files = [files for files in os.listdir(filepath) if 'T' in files]
        Path(output_dir0).mkdir(parents=True, exist_ok=True)
        for file in files:
            output_dir1 = os.path.join(output_dir0, file)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            atoms = read(os.path.join(filepath, file, '%s.traj' % file), '0')
            EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
            try:
                EFzeolite = ExtraFrameworkMaker()
                my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(atoms, EF_atoms)
                write(output_dir1 + '/CuOCu.traj', my_atoms)
            except:
                print(self.sample_zeolite, file, ' failed!')


if __name__ == '__main__':
    topo = ['CHA', 'MOR', 'MWW', 'FAU', 'FER', 'MFI', 'BEC', 'MON', 'MSE', 'AFO', 'AHT', 'BOG', 'CFI',
            'CGF', 'DON', 'EMT', 'EON', 'EUO', 'GON', 'IWS', 'LTF', 'MTW', 'OBW', 'OSI', 'RHO', 'RSN',
            'SBE', 'SSY', 'TER', 'VFI', 'WEI', 'ABW', 'ACO', 'AET', 'AFI', 'AFN', 'AFR', 'AFT', 'AFY',
            'ANA', 'APC', 'APD', 'AST', 'ASV', 'ATN', 'ATO', 'ATT', 'ATV', 'AWW', 'BCT', 'BIK', 'BOF',
            'BRE', 'BSV', 'CAN', 'CGS', 'CHI', 'CLO', 'CON', 'CZP', 'DDR', 'DFO', 'DFT', 'EAB', 'EDI',
            'ERI', 'ESV', 'ETR', 'EZT', 'FAR', 'FRA', 'GIS', 'GIU', 'GME', 'GOO', 'HEU', 'IFR', 'IHW',
            'ISV', 'ITE', 'ITH', 'ITR', 'IWR', 'IWV', 'JBW', 'KFI', 'LEV', 'LIO', 'LOS', 'LOV', 'LTA',
            'LTL', 'LTN', 'MAR', 'MAZ', 'MEL', 'MEP', 'MER', 'MFS', 'MOZ', 'MTF', 'MTN', 'NAB', 'NAT',
            'NES', 'NON', 'NPO', 'NSI', 'OFF', 'OSO', 'PAU', 'PHI', 'PON', 'RRO', 'RTH', 'RUT', 'RWR',
            'RWY', 'SAO', 'SAS', 'SAV', 'SBN', 'SBS', 'SBT', 'SFE', 'SFF', 'SFG', 'SFH', 'SFN', 'SFO',
            'SFS', 'SIV', 'SOD', 'STF', 'STW', 'THO', 'TOL', 'USI', 'UTL', 'VET', 'VNI', 'VSV', 'WEN',
            'YUG', 'ZON', 'IWW', 'STT', 'SVR', 'TUN', 'STO', 'AEI', 'AEL', 'AEN', 'AFG', 'AFS', 'AFX',
            'ATS', 'AVL', 'AWO', 'BOZ', 'BPH', 'CAS', 'CDO', 'DAC', 'DOH', 'EPI', 'FER', 'UWY', 'TON',
            'TSC', 'UEI', 'UFI', 'UOS', 'UOZ', 'SZR', 'STI', 'SVV', 'SGT', 'SOF', 'SOS', 'SSF', 'SAT',
            'SAF', 'RTE', 'PUN', 'PCR', 'OWE', 'PAR', 'NPT', 'MVY', 'MSO', 'MEI', 'LIT', 'LAU', 'LTJ',
            'JOZ', 'JRY', 'JSN', 'JST', 'JSW', 'ITW', 'ITV', 'IRR', 'IMF', 'BEA', 'SVY']
    """
    for sample_zeolite in topo:
        traj_saver = SaveExtraFrameworkConfigurations(sample_zeolite)
        traj_saver.save_bare_zeo(output_dir='/Users/jiaweiguo/Box/00_bare_zeo/')
        traj_saver.save_all_ZCu(dir_name='/Users/jiaweiguo/Box/01_1Al_zeo'):
        print(sample_zeolite + 'is done!')
        
    for sample_zeolite in topo:
        try:
            zeo_dir = '/Users/jiaweiguo/Box/00_bare_zeo/00_%s/%s.traj' % (sample_zeolite, sample_zeolite)
            output_1Al = '/Users/jiaweiguo/Box/all_zeo_database/1Al/00_%s' % sample_zeolite
            output_2Al = '/Users/jiaweiguo/Box/all_zeo_database/2Al/00_%s' % sample_zeolite
            traj_saver = SaveExtraFrameworkConfigurations(sample_zeolite)
            traj_saver.save_all_Al_zeo(zeo_dir, output_1Al, output_2Al)
            print(sample_zeolite + ' is done!')
        except:
            print(sample_zeolite + ' is failed!')
    """
    for sample_zeolite in topo:
        filepath = '/Users/jiaweiguo/Box/all_zeo_database/2Al/00_%s/' % sample_zeolite
        output_dir = '/Users/jiaweiguo/Box/all_zeo_database/2Al_CuOCu/01_%s/' % sample_zeolite
        traj_saver = SaveExtraFrameworkConfigurations(sample_zeolite)
        traj_saver.save_all_CuOCu(filepath, output_dir)
        print(sample_zeolite, ' is done!')
        
