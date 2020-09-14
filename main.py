from source.zeotype import Zeotype
from ase.visualize import view
import os
if __name__ == '__main__':
    from ase.io import read
    cif_path = os.path.join(os.getcwd(), 'source', 'BEA.cif')
    b = read(cif_path)
    z = Zeotype(b)
    index = z.add_cluster(187, 5)
    view(z)
    view(z.clusters[index])

