from source.zeotype import Zeotype
if __name__ == '__main__':
    from ase.io import read
    b = read('.source/BEA.cif')
    z = Zeotype(b)
    z.add_cluster(35, 30)
