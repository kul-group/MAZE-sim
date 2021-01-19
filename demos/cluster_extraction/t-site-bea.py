
from maze import Zeotype, OpenDefect
from ase.visualize import view
from collections import defaultdict

def main():
    cif_dir = "/Users/dda/Code/zeotype/data/BEA.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    open_defect = OpenDefect(zeolite)

    unique_t_site_indices = {}
    for key, value in zeolite.site_to_atom_indices.items():
        if 'T' in key:
            unique_t_site_indices[key] = value[1]
    unique_t_site_to_od = defaultdict(list)
    for key, t_site in unique_t_site_indices.items():

        for o_index in open_defect.neighbor_list.get_neighbors(t_site)[0]:
            for si_index in open_defect.neighbor_list.get_neighbors(o_index)[0]:
                if si_index == t_site:
                    continue
                new_od = open_defect.delete_atoms([si_index]).cap_atoms()
                unique_t_site_to_od[key].append(new_od)
                new_t_site_index = new_od.index_mapper.get_index(open_defect.name, new_od.name, t_site)
                new_od[new_t_site_index].symbol = 'Co'
                view(new_od)
        break


if __name__ == "__main__":
    main()