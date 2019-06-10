#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 16:36:06 2019

@author: ritwit
"""

from pygenarris import *
import numpy as np

lattice = np.eye(3)
xc = np.array([1, 1.3, 1.6]);
yc = np.array([1, 1.3, 1.6]);
zc = np.array([1, 1.3, 1.6]);
atoms = "C C C "
xtal = crystal();
vdw_matrix = np.array([[0, 1, 1],[1 , 0, 1],[ 1, 1, 0]], dtype = 'float32')


create_crystal_from_array(xtal, lattice, xc, yc, zc, atoms,3, 3, 1 )
check_structure_with_vdw_matrix(xtal,  vdw_matrix)