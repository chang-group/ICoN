#!/usr/bin/env python
# coding: utf-8

import pytraj as pt
import numpy as np

TOP = '../../../AB13.prmtop'
TRAJ = './aB13SyntheticE.dcd'
traj = pt.iterload(TRAJ, TOP, stride=1)
print(traj.shape)
