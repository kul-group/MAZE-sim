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


def save_zeolites(sample_zeolite: str):
    EFzeolite = ExtraFrameworkMaker(iza_code=sample_zeolite)

    # 1Al and 2Al replacement
    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    print(EFzeolite.t_site_indices)
    print(EFzeolite.t_site_indices_count)
    print('Number of unique Al pairs: ', len(EFzeolite.traj_2Al))

    # sample all Bronsted sites
    dict_1Al_ZH = EFzeolite.get_all_Bronsted_sites(case_1Al=True)
    print('Single Al Bronsted site sampling is done!')
    dict_2Al_ZH = EFzeolite.get_all_Bronsted_sites(case_2Al=True)
    print('Two Al Bronsted site sampling is done!')

    output_dir = '/Users/jiaweiguo/Box/all_zeo_database'

    for site_name, traj in dict_1Al_ZH.items():
        output_dir1 = os.path.join(output_dir, sample_zeolite, '1Al', site_name, 'H')
        Path(output_dir1).mkdir(parents=True, exist_ok=True)

        for count, atoms in enumerate(traj):
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = nl.get_neighbors(index_H[0])[0]
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val != index_H][0]
            output_filename = sample_zeolite + '_' + '1Al_' + site_name + '_H_O' + str(neigh_o_indices)
            # use O tag instead of indices (maybe)
            my_path = os.path.join(output_dir1, output_filename)
            write(my_path+'.traj', atoms)

    for site_name, traj in dict_2Al_ZH.items():
        output_dir2 = os.path.join(output_dir, sample_zeolite, '2Al', site_name, 'H')
        Path(output_dir2).mkdir(parents=True, exist_ok=True)

        for atoms in traj:
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = list(nl.get_neighbors(index_H[0])[0]) + list(nl.get_neighbors(index_H[1])[0])
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val not in index_H]
            assert len(neigh_o_indices) == 2
            output_filename = sample_zeolite + '_' + '2Al_H_' + site_name + '_O' + str(neigh_o_indices[0]) + '_O' + str(neigh_o_indices[1])
            my_path = os.path.join(output_dir2, output_filename)
            # write(my_path, atoms, format='vasp') write as poscar
            write(my_path+'.traj', atoms)

    # insert extra-framework CuOCu
    EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
    EF_type = 'CuOCu'
    for site_name, zeolite in EFzeolite.dict_2Al_replaced.items():
        output_dir3 = os.path.join(output_dir, sample_zeolite, '2Al', site_name, EF_type)
        Path(output_dir3).mkdir(parents=True, exist_ok=True)

        atoms = EFzeolite.insert_ExtraFrameworkAtoms(zeolite, EF_atoms)
        output_filename = sample_zeolite + '_' + '2Al_CuOCu_' + site_name
        my_path = os.path.join(output_dir3, output_filename)
        write(my_path+'.traj', atoms)

    print('Extra-framework atoms insertion is done!')

    print('File writing is finished!')
    # ase.build('CuOCu')
    #TODO: ADD TAGS FOR O AS WELL (BESIDES INDEX)
    #TODO: ADD FOLDER STRUCTURE FUNCTION IN


def save_zeolites_v2(sample_zeolite: str, output_dir):
    # save bare zeolite structure
    EFzeolite = ExtraFrameworkMaker(iza_code=sample_zeolite)
    output_dir0 = os.path.join(output_dir, sample_zeolite)
    Path(output_dir0).mkdir(parents=True, exist_ok=True)
    my_path = os.path.join(output_dir0, sample_zeolite)
    write(my_path + '.traj', EFzeolite.EFzeolite)

    # 1Al and 2Al replacement
    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    print(EFzeolite.t_site_indices)
    print(EFzeolite.t_site_indices_count)
    print('Number of unique Al pairs: ', len(EFzeolite.traj_2Al))

    for site_name, atoms_1Al in EFzeolite.dict_1Al_replaced.items():
        output_dir1 = os.path.join(output_dir, sample_zeolite, '1Al', site_name)
        Path(output_dir1).mkdir(parents=True, exist_ok=True)
        my_path = os.path.join(output_dir1, sample_zeolite + '_' + '1Al_' + site_name)
        write(my_path + '.traj', atoms_1Al[0])  # only save one structure for each unique T

    for site_name, atoms_2Al in EFzeolite.dict_2Al_replaced.items():
        output_dir2 = os.path.join(output_dir, sample_zeolite, '2Al', site_name)
        Path(output_dir2).mkdir(parents=True, exist_ok=True)
        my_path = os.path.join(output_dir2, sample_zeolite + '_' + '2Al_' + site_name)
        write(my_path + '.traj', atoms_2Al)

    """
    name_out = 'Hb_%s_%s_%s_%s' % (topo,str(i_tag).zfill(2),str(tag_T).zfill(2),str(count).zfill(2))
        #view(atoms_new)
        
        c = FixAtoms(indices=[atom.index for atom in atoms_new if atom.symbol != 'H'])
        atoms_new.set_constraint(c)
        
        if not os.path.exists(new_dir+'/'+name_out):
            os.mkdir(new_dir+'/'+name_out)
        io.write(new_dir+'/'+name_out+'/'+name_out+'.traj',atoms_new)
        io.write(new_dir+'/'+name_out+'/'+'manual_start'+'.traj',atoms_new)
        print(name_out)
"""


def save_H_structures(sample_zeolite, Al_num: str):
    folder = '/Users/jiaweiguo/Box/all_zeo_database/' + sample_zeolite + '/' + Al_num
    filepaths = [os.path.join(folder, dirs) for dirs in os.listdir(folder) if 'T' in dirs]
    for filepath in filepaths:
        files = [files for files in os.listdir(filepath) if '.traj' in files]
        for file in files:
            atoms = read(os.path.join(filepath, file))
            EFzeolite = ExtraFrameworkMaker()
            dict_H = EFzeolite.get_Bronsted_sites(atoms)
            for site_name, my_atoms in dict_H.items():
                output_dir1 = os.path.join(filepath, 'H')
                Path(output_dir1).mkdir(parents=True, exist_ok=True)
                my_path = os.path.join(output_dir1, sample_zeolite + '_' + Al_num + '_' + site_name)
                write(my_path + '.traj', my_atoms)


def save_EF_structures(sample_zeolite, EF_atoms, Al_num: str, EF_tag: str):
    folder = '/Users/jiaweiguo/Box/all_zeo_database/' + sample_zeolite + '/' + Al_num
    filepaths = [os.path.join(folder, dirs) for dirs in os.listdir(folder) if 'T' in dirs]
    for filepath in filepaths:
        files = [files for files in os.listdir(filepath) if '.traj' in files]
        for file in files:
            zeolite = read(os.path.join(filepath, file))
            EFzeolite = ExtraFrameworkMaker()
            atoms = EFzeolite.insert_ExtraFrameworkAtoms(zeolite, EF_atoms)
            output_dir1 = os.path.join(filepath, EF_tag)
            Path(output_dir1).mkdir(parents=True, exist_ok=True)
            my_path = os.path.join(output_dir1, sample_zeolite + '_' + Al_num + '_' + EF_tag)
            write(my_path + '.traj', atoms)


def get_all(zeo_code):
    # save_zeolites('CHA') older version
    save_zeolites_v2(zeo_code, '/Users/jiaweiguo/Box/all_zeo_database')
    save_H_structures(zeo_code, '1Al')
    save_H_structures(zeo_code, '2Al')
    EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
    EF_type = 'CuOCu'
    save_EF_structures(zeo_code, EF_atoms, '2Al', EF_type)


def save_bare_zeo(sample_zeolite):
    output_dir = '/Users/jiaweiguo/Box/00_bare_zeo/'
    EFzeolite = ExtraFrameworkMaker(iza_code=sample_zeolite)
    output_dir0 = os.path.join(output_dir, '00_'+sample_zeolite)
    Path(output_dir0).mkdir(parents=True, exist_ok=True)
    my_path = os.path.join(output_dir0, sample_zeolite)
    write(my_path + '.traj', EFzeolite.EFzeolite)


def save_all_Al_zeo(sample_zeolite):
    EFzeolite = ExtraFrameworkMaker(sample_zeolite, '/Users/jiaweiguo/Box/00_test_zeo/00_%s/new_%s.traj'
                                    % (sample_zeolite, sample_zeolite))

    # 1Al and 2Al replacement
    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    print(EFzeolite.t_site_indices)
    print(EFzeolite.t_site_indices_count)
    print('Number of unique Al pairs: ', len(EFzeolite.traj_2Al))

    output_dir = '/Users/jiaweiguo/Box/01_1Al_zeo/01_%s' % sample_zeolite
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for site_name, atoms_1Al in EFzeolite.dict_1Al_replaced.items():
        output_dir1 = os.path.join(output_dir, site_name)
        Path(output_dir1).mkdir(parents=True, exist_ok=True)
        my_path = os.path.join(output_dir1, site_name)
        write(my_path + '.traj', atoms_1Al[0])  # only save one structure for each unique T

    output_dir = '/Users/jiaweiguo/Box/01_2Al_zeo/01_%s' % sample_zeolite
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for site_name, atoms_2Al in EFzeolite.dict_2Al_replaced.items():
        output_dir2 = os.path.join(output_dir, site_name)
        Path(output_dir2).mkdir(parents=True, exist_ok=True)
        my_path = os.path.join(output_dir2, site_name)
        write(my_path + '.traj', atoms_2Al)


def save_all_Z_TM(sample_zeolite):
    original_dir = '/Users/jiaweiguo/Box/01_1Al_zeo'
    for dir_name in os.listdir(original_dir):
        if sample_zeolite in dir_name:
            working_dir = original_dir + '/' + dir_name
            os.chdir(working_dir)
            for sub_dir_name in os.listdir(working_dir):
                if 'T' in sub_dir_name:
                    sub_working_dir = working_dir + '/' + sub_dir_name
                    os.chdir(sub_working_dir)
                    track_list = []
                    for root, dirs, files in os.walk(sub_working_dir):
                        [track_list.append(dir_names) for dir_names in dirs if "opt_" in dir_names]
                        break  # only track one layer down the folder
                    if len(track_list) != 0:
                        track_list.sort(reverse=True)
                        atoms = read(track_list[0] + "/vasprun.xml", '-1')
                        EFzeolite = ExtraFrameworkMaker()
                        d_Z_Cu = 2.6
                        dict_ZCu = EFzeolite.get_Z_TM(atoms, d_Z_Cu, 'Cu')

                        output_dir = '/Users/jiaweiguo/Box/02_1Cu_zeo/01_%s' % sample_zeolite
                        if not os.path.exists(output_dir):
                            os.mkdir(output_dir)

                        for site_name, my_atoms in dict_ZCu.items():
                            output_dir1 = os.path.join(output_dir, sub_dir_name + '_' + site_name)
                            Path(output_dir1).mkdir(parents=True, exist_ok=True)
                            nl = NeighborList(natural_cutoffs(my_atoms), bothways=True, self_interaction=False)
                            nl.update(my_atoms)
                            c = FixAtoms(indices=[atom.index for atom in my_atoms if atom.symbol not in ['Cu']])
                            my_atoms.set_constraint(c)
                            write(output_dir1 + '/1Cu.traj', my_atoms)

                    os.chdir(sub_working_dir)

            os.chdir(original_dir)


def save_all_CuOCu(sample_zeolite):
    original_dir = '/Users/jiaweiguo/Box/01_2Al_zeo'
    for dir_name in os.listdir(original_dir):
        if sample_zeolite in dir_name:
            working_dir = original_dir + '/' + dir_name
            os.chdir(working_dir)
            for sub_dir_name in os.listdir(working_dir):
                if 'T' in sub_dir_name:
                    sub_working_dir = working_dir + '/' + sub_dir_name
                    os.chdir(sub_working_dir)
                    track_list = []
                    for root, dirs, files in os.walk(sub_working_dir):
                        [track_list.append(dir_names) for dir_names in dirs if "opt_" in dir_names]
                        break  # only track one layer down the folder
                    if len(track_list) != 0:
                        track_list.sort(reverse=True)
                        atoms = read(track_list[0] + "/vasprun.xml", '-1')

                        EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
                        EF_type = 'CuOCu'
                        EFzeolite = ExtraFrameworkMaker()
                        my_atoms = EFzeolite.insert_ExtraFrameworkAtoms(atoms, EF_atoms)

                        output_dir = '/Users/jiaweiguo/Box/02_2Cu_zeo/01_%s' % sample_zeolite
                        if not os.path.exists(output_dir):
                            os.mkdir(output_dir)

                        output_dir1 = os.path.join(output_dir, sub_dir_name)
                        Path(output_dir1).mkdir(parents=True, exist_ok=True)
                        nl = NeighborList(natural_cutoffs(my_atoms), bothways=True, self_interaction=False)
                        nl.update(my_atoms)

                        Index = [a.index for a in my_atoms if a.symbol in ['Al', 'Cu']]
                        neigh_list = []
                        for index in Index:
                            indices, offsets = nl.get_neighbors(index)
                            neigh_list.extend(indices)
                        print(neigh_list)
                        c = FixAtoms(indices=[atom.index for atom in my_atoms if atom.symbol not in ['Cu', 'Al'] and
                                              atom.index not in neigh_list])
                        my_atoms.set_constraint(c)
                        write(output_dir1 + '/2Cu.traj', my_atoms)

                    os.chdir(sub_working_dir)

            os.chdir(original_dir)


def make_bare_zeo():
    """
    zeo = ['CHA', 'SOD', 'LTA', 'FAU', 'MOR', 'MFI', 'MAZ', 'BEA', 'AEI', 'RHO', 'MWW', 'SVY']
    for each_zeo in zeo:
        get_all(each_zeo)
        print(each_zeo + 'is done!')
    """
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
            'JOZ', 'JRY', 'JSN', 'JST', 'JSW', 'ITW', 'ITV', 'IRR', 'IMF']
    print(len(topo))

    addition = ['BEA', 'SVY']
    for each_zeo in addition:
        save_bare_zeo(each_zeo)


if __name__ == '__main__':
    """
    zeo = ['CHA', 'SOD', 'LTA', 'FAU', 'MOR', 'MFI', 'MAZ', 'BEA', 'AEI', 'RHO', 'MWW']
    for each_zeo in zeo:
        # save_all_Z_TM(each_zeo)
        # save_all_Al_zeo(each_zeo)
        save_all_CuOCu(each_zeo)
        print(each_zeo + ' is done!')
    """
    save_all_CuOCu('CHA')
