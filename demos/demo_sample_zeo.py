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
from xlsxwriter import Workbook


def sample_all_zeo_with_certain_T_sites(topo, t_count):
    output = []
    for sample_zeolite in topo:
        EFzeolite = ExtraFrameworkMaker(iza_code=sample_zeolite)
        EFzeolite.get_t_sites()
        if len(EFzeolite.t_site_indices) == t_count:
            output.append([sample_zeolite, len(EFzeolite.EFzeolite)])
    return output


if __name__ == '__main__':
    
    topo = ['CHA', 'MOR', 'MWW', 'FAU', 'SOD', 'MFI', 'RHO', 'AEI', 'LTA', 'MAZ', 'BEA', 'BOG', 'CFI', 
            'CGF', 'DON', 'EMT', 'EON', 'EUO', 'GON', 'IWS', 'LTF', 'MTW', 'OBW', 'OSI', 'BEC', 'RSN',
            'SBE', 'SSY', 'TER', 'VFI', 'WEI', 'ABW', 'ACO', 'AET', 'AFI', 'AFN', 'AFR', 'AFT', 'AFY',
            'ANA', 'APC', 'APD', 'AST', 'ASV', 'ATN', 'ATO', 'ATT', 'ATV', 'AWW', 'BCT', 'BIK', 'BOF',
            'BRE', 'BSV', 'CAN', 'CGS', 'CHI', 'CLO', 'CON', 'CZP', 'DDR', 'DFO', 'DFT', 'EAB', 'EDI',
            'ERI', 'ESV', 'ETR', 'EZT', 'FAR', 'FRA', 'GIS', 'GIU', 'GME', 'GOO', 'HEU', 'IFR', 'IHW',
            'ISV', 'ITE', 'ITH', 'ITR', 'IWR', 'IWV', 'JBW', 'KFI', 'LEV', 'LIO', 'LOS', 'LOV', 'MSE',
            'LTL', 'LTN', 'MAR', 'AFO', 'MEL', 'MEP', 'MER', 'MFS', 'MOZ', 'MTF', 'MTN', 'NAB', 'NAT',
            'NES', 'NON', 'NPO', 'NSI', 'OFF', 'OSO', 'PAU', 'PHI', 'PON', 'RRO', 'RTH', 'RUT', 'RWR',
            'RWY', 'SAO', 'SAS', 'SAV', 'SBN', 'SBS', 'SBT', 'SFE', 'SFF', 'SFG', 'SFH', 'SFN', 'SFO',
            'SFS', 'SIV', 'FER', 'STF', 'STW', 'THO', 'TOL', 'USI', 'UTL', 'VET', 'VNI', 'VSV', 'WEN',
            'YUG', 'ZON', 'IWW', 'STT', 'SVR', 'TUN', 'STO', 'MON', 'AEL', 'AEN', 'AFG', 'AFS', 'AFX',
            'ATS', 'AVL', 'AWO', 'BOZ', 'BPH', 'CAS', 'CDO', 'DAC', 'DOH', 'EPI', 'FER', 'UWY', 'TON',
            'TSC', 'UEI', 'UFI', 'UOS', 'UOZ', 'SZR', 'STI', 'SVV', 'SGT', 'SOF', 'SOS', 'SSF', 'SAT',
            'SAF', 'RTE', 'PUN', 'PCR', 'OWE', 'PAR', 'NPT', 'MVY', 'MSO', 'MEI', 'LIT', 'LAU', 'LTJ',
            'JOZ', 'JRY', 'JSN', 'JST', 'JSW', 'ITW', 'ITV', 'IRR', 'IMF', 'AHT']
    
    # print(sample_all_zeo_with_certain_T_sites(topo, t_count=1))

    wb = Workbook('/Users/jiaweiguo/Box/temp_T_count.xlsx')
    sheet1 = wb.add_worksheet('Sheet1')

    output = []
    for count, zeolite in enumerate(topo):
        try:
            EFzeolite = ExtraFrameworkMaker(iza_code=zeolite)
            EFzeolite.get_t_sites()
            output.append([zeolite, len(EFzeolite.t_site_indices), len(EFzeolite.EFzeolite)])
            sheet1.write(count, 0, zeolite)
            sheet1.write(count, 1, len(EFzeolite.t_site_indices))
            sheet1.write(count, 2, len(EFzeolite.EFzeolite))
            print([zeolite, len(EFzeolite.t_site_indices), len(EFzeolite.EFzeolite)])
        except:
            print([zeolite, 0, 0])

    wb.close()
