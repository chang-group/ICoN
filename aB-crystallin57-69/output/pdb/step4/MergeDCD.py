#!/usr/bin/env python

import pytraj as pt
import sys
import os

TOP = '../../../top.prmtop'
path = './'
DIR_SAVE = './'
file_name = 'AB13'

ls_files = os.listdir(path) # your directory path
n_files = len(ls_files)     # n files 
TRAJ = [path + i for i in ls_files]

print(TRAJ)

traj = pt.iterload(TRAJ, TOP, stride=1)
print(traj.shape)

pt.write_traj(DIR_SAVE + file_name + 'SyntheticNovel' + '.dcd', traj, overwrite=True)
