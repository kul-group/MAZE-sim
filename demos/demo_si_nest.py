from maze import Zeotype
from ase.visualize import view


def main():
    cif_dir = "//data/GOO.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    view(zeolite)
    iz = zeolite.get_imperfect_zeotype()
    iz = iz.create_silanol_defect(95)
    view(iz)
    # view(open_framework)
    # open_framework.remove_caps()
    # view(open_framework)
    # open_framework.integrate_other_zeotype(cluster)
    # view(open_framework)


if __name__ == '__main__':
    main()
