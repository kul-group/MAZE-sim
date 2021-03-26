from maze.zeolite import Zeolite
from ase.visualize import view


def main():
    cha = Zeolite.make('CHA')
    cluster, od = cha.get_cluster(10)
    # print(cluster.additions)
    # cluster = cluster.remove_caps('h_caps', 'h_caps_6_')
    # print(cluster.name)
    # print(cluster.index_mapper.main_index)
    # open_framework = open_framework.cap_atoms()
    # #view(open_framework)
    # view(cluster)
    #
    # pass
    # #cluster = cluster.cap_atoms()

    # view(cluster)
    # view(open_framework.cap_atoms())
    # view(open_framework)

    #open_framework = open_framework.remove_caps()
    # view(open_framework)
    # open_framework = open_framework.integrate_other_zeotype(cluster)
    # view(open_framework)


if __name__ == '__main__':
    main()
