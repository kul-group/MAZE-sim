from source.zeotype import Zeotype, ImperfectZeotype
from ase.visualize import view


def main():
    cif_dir = "/Users/dda/Code/zeotype/data/GOO.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    view(zeolite)
    iz = zeolite.get_imperfect_zeolite()
    iz = iz.create_silanol_defect(95)
    view(iz)
    # view(open_framework)
    # open_framework.remove_caps()
    # view(open_framework)
    # open_framework.integrate_other_zeotype(cluster)
    # view(open_framework)


if __name__ == '__main__':
    main()
