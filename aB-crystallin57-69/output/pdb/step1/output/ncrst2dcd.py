#!/usr/bin/env python

import pytraj as pt
import sys
import os

TOP = '../../top.prmtop'
path = './ncrst/'
DIR_SAVE = './'
file_name = '1z0q_1SyntheticAll' # file name
ls_files = os.listdir(path) # your directory path
n_files = len(ls_files)     # n files 

TRAJ = [path + 'Min_' + str(i) + '.ncrst' for i in range(1,n_files+1)]
traj = pt.iterload(TRAJ, TOP, stride=1)
pt.write_traj(DIR_SAVE + file_name + '.dcd', traj, overwrite=True)
