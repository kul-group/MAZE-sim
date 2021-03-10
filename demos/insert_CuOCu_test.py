from maze.zeotypes import Zeotype, ImperfectZeotype
from maze.extra_framework import ExtraFramework
from collections import defaultdict
from ase.neighborlist import natural_cutoffs, NeighborList
from ase import Atoms
import numpy as np
from ase.visualize import view
import copy as copy
from ase.io import write, read
import ase.db

# MFI zeolite
db = ase.db.connect('/Users/jiaweiguo/Box/CuOCu_systematic/systematic_cuocu.only_final.db')
e_zcu = {1: -2303.65016,
         2: -2303.55322,
         3: -2303.711655,
         4: -2303.950027,
         5: -2303.851314,
         6: -2303.59494,
         7: -2303.685289,
         8: -2303.715035,
         9: -2303.632148,
         10: -2303.669027,
         11: -2303.890084,
         12: -2303.736395, }
e_o2 = -9.88984331
e_sili_MFI = -2303.68129361

traj = []
Al_pair = []
repeated_pair = []
E_repeated = []
d_AlAl_list = []
for row in db.select(min_energy_structure=True, tag_al1=1, tag_al2=5):
    Al1_tag = row.tag_al1
    Al2_tag = row.tag_al2
    atoms = db.get_atoms(id=row.id)
    traj.append(atoms)
    e = row.energy + e_sili_MFI - e_zcu[Al1_tag] - e_zcu[Al2_tag] - 0.5 * e_o2
    repeated_pair = [[Al1_tag, Al2_tag]]
    E_repeated = [e]

    index_Al = [a.index for a in atoms if a.symbol == 'Al']
    d_AlAl = atoms.get_distance(index_Al[0], index_Al[1])
    d_AlAl_list.append(d_AlAl)
    if [Al1_tag, Al2_tag] not in Al_pair and [Al2_tag, Al1_tag] not in Al_pair:
        Al_pair.append([Al1_tag, Al2_tag])
    else:
        repeated_pair.append([Al1_tag, Al2_tag])
        E_repeated.append(e)

print(len(traj))
# print(Al_pair)
# print(len(Al_pair))
print(repeated_pair)
# view(traj)
print(E_repeated)
print(d_AlAl_list)
'''
cif_dir = 'MFI.cif'
my_zeolite = ExtraFramework(cif_dir=cif_dir)
# unique_t_site_indices = my_zeolite.get_unique_t_sites(cif_dir)
my_zeolite.create_1Al_replacement()
my_zeolite.create_2Al_replacement()
# write('BEA_2Al_replaced.traj', my_zeolite.traj_2Al)
print(my_zeolite.count)
'''
## 88 vs 69 (remove repeated Al-Al pair)
