#!/usr/bin/env python
# coding: utf-8

import pytraj as pt
import numpy as np
#import pandas as pd

TOP = '../../../top.prmtop'
TRAJ = './1z0q_1SyntheticE.dcd'
traj = pt.iterload(TRAJ, TOP, stride=1)
print(traj.shape)
