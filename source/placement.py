from sklearn.cluster import MeanShift
from scipy.optimize import minimize
from ase import Atom
import numpy as np
from random import random

# TODO: reduce number of input parameters -> more automatic

def min_dist(pos, host):
# returns the distance to the closest atom in the host at position
    dummy_atom = Atom('H', position=pos)
    dummy_host = host + dummy_atom
    return(min(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True)))

def avg_dist(pos, host):
# returns the average distance to atoms in the host
    dummy_atom = Atom('H', position=pos)
    dummy_host = host + dummy_atom
    return(np.average(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True)))

def find_void(pos, host):
# from starting position, finds the nearest void in the host structure
    guess = pos
    func = lambda pos: -1*min_dist(pos, host) # 1 param function for scipy.minimize
    ans = minimize(func, guess)
    return(ans.x)

def sphere_sample(radius, num_pts=None):
# generates random positions on the surface of a sphere of certain radius
    if num_pts == None:
        num_pts = 300
    vect_list = []
    for i in range(num_pts):
        x, y, z = [2*random()-1 for i in range(3)]
        if x**2 + y**2 + z**2 > 1.0:
            pass
        else:
            unit_vec = [x, y, z]/np.linalg.norm([x,y,z])
            vect_list.append(radius*unit_vec)
    return(vect_list)

def get_place_clusters(host, index, radius, num_pts, cutoff):
# rejects random points around host atom which are too close to another host atom
# clusters the non-rejected points
# if num_pts is too small, will get error
    assert (radius > cutoff)
    guess_pos = sphere_sample(radius, num_pts)
    viable_pos = []
    host_pos = host.get_positions()[index]
    for i in guess_pos:
        d = min_dist(i+host_pos, host) # cutoff must be smaller than radius
        if d > cutoff:
            viable_pos.append(i+host_pos)
    ms = MeanShift(bin_seeding=True)
    ms.fit(np.array(viable_pos))
    cluster_centers = ms.cluster_centers_
    return(cluster_centers)

def find_best_place(host, index, radius, num_pts, cutoff):
# picks the best location to place an adsorbate around the host atom

    viable_pos = get_place_clusters(host, index, radius, num_pts, cutoff)
    best_avg_dist = 0
    for pos in viable_pos:
        void = find_void(pos, host)
        void_avg_dist = avg_dist(void, host)  # average dist at void nearest to pos
        if void_avg_dist > best_avg_dist:
            best_avg_dist = void_avg_dist     # selects pos with largest nearby void
            best_pos = pos
    return(best_pos)

# testing
if __name__=='__main__':

    from ase.io import read
    from ase.visualize import view
    host = read('BEA.cif')
    a = get_place_clusters(host, 185, 2.2, 400, 1.9)
    viz = host.copy()
    c = []
    best_pos = find_best_place(host, 185, 2.2, 400, 1.9)
    print(best_pos)
    viz = viz + Atom('H', position=best_pos) # visualization
    view(viz)