#!/usr/bin/env python

import pytraj as pt
import sys
import os

TOP = '../../top.prmtop'
path = './ncrst/'
DIR_SAVE = './'
file_name = '2nao_4SyntheticAll' # input file name,
                                 # all DL generated synthetic confs

ls_files = os.listdir(path) # path to minimized .ncrst confs
n_files = len(ls_files)     # number of files 

TRAJ = [path + 'Min_' + str(i) + '.ncrst' for i in range(1,n_files+1)]
traj = pt.iterload(TRAJ, TOP, stride=1)
pt.write_traj(DIR_SAVE + file_name + '.dcd', traj, overwrite=True)
