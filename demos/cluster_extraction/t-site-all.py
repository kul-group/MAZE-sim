from ase.io import  write
from maze import Zeotype, OpenDefect
import os
from pathlib import Path
from collections import defaultdict
from maze import download_cif

def defect_maker(zeolite_code, output_dir):
    cif_folder = os.path.join(output_dir, 'cif_files')
    download_cif(zeolite_code, data_dir=cif_folder)
    cif_dir = os.path.join(cif_folder, zeolite_code + '.cif')

    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    open_defect = OpenDefect(zeolite)
    unique_t_site_indices = {}
    for site_name, value in zeolite.site_to_atom_indices.items():
        if 'T' in site_name:
            unique_t_site_indices[site_name] = value[1]

    # build dictionary of open defects
    unique_t_site_to_od = defaultdict(list)
    for site_name, t_site in unique_t_site_indices.items():

        for o_index in open_defect.neighbor_list.get_neighbors(t_site)[0]:
            for si_index in open_defect.neighbor_list.get_neighbors(o_index)[0]:
                if si_index == t_site:
                    continue
                new_od = open_defect.delete_atoms([si_index]).cap_atoms()
                unique_t_site_to_od[site_name].append(new_od)


    # save T sites
    for site_name, od_list in unique_t_site_to_od.items():
        output_dir = os.path.join(output_dir, zeolite_code, site_name)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        for index, od in enumerate(od_list):
            output_filename = zeolite_code + '_' + site_name + '_' + str(index) + '.traj'
            my_path = os.path.join(output_dir, output_filename)
            write(my_path, od)


if __name__ == "__main__":
    defect_maker('BEA', '/Users/dda/Desktop/site')
