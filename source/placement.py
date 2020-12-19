from sklearn.cluster import MeanShift
from scipy.optimize import minimize
from ase import Atom
import numpy as np
from random import random

def min_dist(pos, host):
    """minimum distance from position to a host atom :param pos: vector x,y,z
    position :param host: host atoms object :return: float, minimum distance

    Args:
        pos:
        host:
    """
    dummy_atom = Atom('H', position=pos)
    dummy_host = host + dummy_atom
    return(min(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True)))

def avg_dist(pos, host):
    """average distance from position to all host atoms :param pos: vector x,y,z
    position :param host: host atoms object :return: float, average distance

    Args:
        pos:
        host:
    """
    dummy_atom = Atom('H', position=pos)
    dummy_host = host + dummy_atom
    return(np.average(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True)))

def find_void(pos, host):
    """finds nearest empty region in host :param pos: vector x,y,z position
    :param host: host atoms object :return: vector position of center of empty
    void

    Args:
        pos:
        host:
    """
    guess = pos
    func = lambda pos: -1*min_dist(pos, host) # 1 param function for scipy.minimize
    ans = minimize(func, guess)
    return(ans.x)

def sphere_sample(radius, num_pts=None):
    """generates random positions on the surface of a sphere of certain radius
    :param radius: radius of sphere surface to sample :param num_pts: number of
    points to try :return: list of x,y,z positions on surface of sphere

    Args:
        radius:
        num_pts:
    """
    if num_pts == None:
        num_pts = 500
    vect_list = []
    for i in range(num_pts):
        x, y, z = [2*random()-1 for i in range(3)]
        # select points inside of unit sphere
        if x**2 + y**2 + z**2 > 1.0:
            pass
        else:
            unit_vec = [x, y, z]/np.linalg.norm([x,y,z])  # normalize vector length to 1
            vect_list.append(radius*unit_vec)  # lengthen vector to desired radius
    return(vect_list)

def get_place_clusters(host, index, radius, cutoff, num_pts=None):
    """finds positions near host atom far enough from other framework atoms,
    clusters these viable positions and returns the centers of these clusters.
    If number of points is too small will return error :param host: host atoms
    object :param index: index of host atom at center :param radius: radius
    around host atom to test points :param cutoff: minimum distance from other
    host atoms allowed for test points :param num_pts: number of points to try
    :return: list. center positions of clusters of points which meet criteria

    Args:
        host:
        index:
        radius:
        cutoff:
        num_pts:
    """
    assert (radius > cutoff)
    if num_pts == None:
        num_pts = 500

    guess_pos = sphere_sample(radius, num_pts)
    host_pos = host.get_positions()[index]
    viable_pos = []

    for pos in guess_pos:
        dist = min_dist(pos+host_pos, host)
        if dist > cutoff:
            viable_pos.append(pos+host_pos)
    ms = MeanShift(bin_seeding=True)
    ms.fit(np.array(viable_pos))
    cluster_centers = ms.cluster_centers_
    return(cluster_centers)

def find_best_place(host, index, radius, cutoff, num_pts=None):
    """picks the best location to place an adsorbate around the host atom :param
    host: host atoms object :param index: index of host atom at center :param
    radius: radius around host atom to test points :param cutoff: minimum
    distance from other host atoms allowed for test points :param num_pts:
    number of points to try :return: vector. x,y,z position of best location

    Args:
        host:
        index:
        radius:
        cutoff:
        num_pts:
    """
    if num_pts == None:
        num_pts = 500
    viable_pos = get_place_clusters(host, index, radius, cutoff, num_pts)
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
    best_pos = find_best_place(host, 185, 2.2, 1.9)
    viz = host + Atom('H', position=best_pos)  # visualization
    view(viz)