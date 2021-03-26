from maze.zeolite import Zeolite, Cluster


if __name__ == "__main__":
    my_iz = Cluster.make('GOO')
    my_iz = my_iz.delete_atoms([0,1,2,3,4,5,6,7])
    print(my_iz.get_type(1))
    print(my_iz.get_site_type())
    print(my_iz.find_type('O1'))