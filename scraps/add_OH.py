#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 10:10:55 2021

@author: jiaweiguo
"""

import os
from pathlib import Path
from ase.io import write, read
import copy
import shutil
from glob import glob
from ase.constraints import FixAtoms
from collections import defaultdict
from ase.neighborlist import natural_cutoffs, NeighborList, mic
from ase import Atoms
import numpy as np
from ase.visualize import view
import ase.db
import pickle
import matplotlib.pyplot as plt
import math
import re

zeolite, TM_type = 'MFI', 'Co'

tol = 1.5 * 2
folder = '/Users/jiaweiguo/Box/P1_pair_site/%s_1Al_%s' %(zeolite, TM_type)
filepaths = [dirs for dirs in os.listdir(folder) if 'T' in dirs]
output_dir0 = '/Users/jiaweiguo/Box/P1_pair_site/%s_1Al_%sOH' %(zeolite, TM_type)
Path(output_dir0).mkdir(parents=True, exist_ok=True)

for file in filepaths:
    try:
        atoms = read(os.path.join(folder, file) + '/opt_400/opt_from_vasp.traj', '0') 
    except:
        atoms = read(os.path.join(folder, file) + '/opt_400/vasprun.xml', '-1') 
         
    output_dir1 = os.path.join(output_dir0, file)
    Path(output_dir1).mkdir(parents=True, exist_ok=True)
    
    Al_index = [atom.index for atom in atoms if atom.symbol == 'Al']
    TM_index = [atom.index for atom in atoms if atom.symbol == TM_type]
    TM_pos = atoms.get_positions()[TM_index]
    
    vec = np.random.normal(size=(3,))
    oh_pos = [TM_pos[0] + vec * 2, TM_pos[0] + vec * 3]
    oh_atoms = Atoms('OH', positions=oh_pos)
    
    oh_cop = np.sum(oh_atoms.positions, 0)/len(oh_atoms)
    distances = mic(oh_cop - oh_atoms.positions, atoms.cell)
    distances = np.linalg.norm(distances, axis=1)
    
    while min(distances) < tol:
        vec = np.random.normal(size=(3,))
        oh_pos = [TM_pos[0] + vec * 2, TM_pos[0] + vec * 3]
        oh_cop = np.sum(oh_atoms.positions, 0)/len(oh_atoms)
        
    
        oh_atoms.translate(oh_pos - oh_cop)
        oh_cop = np.sum(oh_atoms.positions, 0)/len(oh_atoms)

        distances = mic(oh_cop - oh_atoms.positions, atoms.cell)
        distances = np.linalg.norm(distances, axis=1)

    atoms = atoms + Atoms('OH', positions=oh_pos)

    write(output_dir1 + '/starting.traj', atoms)

    
    
    
    
    
    