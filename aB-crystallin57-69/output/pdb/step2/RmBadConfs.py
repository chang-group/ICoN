#!/usr/bin/env python

# Import Libraries
import numpy as np
import pandas as pd
import sys
import os
import pytraj as pt

thresh = -200.0
name = 'AB13'
# load energies
path = '../' + name + '_energies.csv'
ener_p = pd.read_csv(path, sep=",") 
ener_p.columns = ['conf_id','tot', 'bond', 'angle', 'dihedral', 'vdw', 'elec', 'gb', 'vdw_14','elec_14']

print('Here is how energy values look like!')
print(ener_p[['tot', 'bond', 'angle', 'dihedral', 'vdw', 'elec', 'gb', 'vdw_14','elec_14']].head()) # kca


dum_ints = np.arange(ener_p.shape[0])
sel_frames = dum_ints[ener_p['tot'] < thresh]
print('Percenntage of ramaining conformations -- ',sel_frames.shape[0]/dum_ints.shape[0])
print('Number of conf remaining', sel_frames.shape[0])

TRAJ = name + 'Synthetic' + '.dcd'
TOP = '../../../AB13.prmtop'
traj = pt.iterload(TRAJ, TOP, stride=1)


#output trajs
pt.write_traj(name+ 'SyntheticE.dcd', traj[sel_frames], overwrite=True)
#np.save('KeptConfIds.npy',np.array(sel_frames))
