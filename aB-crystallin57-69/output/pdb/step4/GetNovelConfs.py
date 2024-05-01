#!/usr/bin/env python
import pytraj as pt
import numpy as np
import sys
import os


def FindCloseConfs(start, end, traj, traj_s, cut=1.2, mask=['@N,@CA,@C,@O']):

    count = 0
    indx_match = []
    indx_no_match = []

    n_frames2 = traj_s.shape[0]

    for i in range(start, end):

        print('Processing conf -- ', i)
        p_rmsd = pt.rmsd(traj, ref=traj_s[i], mask=mask)

        if (p_rmsd < cut).sum() > 0:
            count += 1
            indx_match.append(i)
            print("Match found at conformation id -- ", i)
        else:
            indx_no_match.append(i)
            print("Novel conf found at conformation id -- ", i)

    return indx_match, indx_no_match, count/n_frames2, count


if __name__ == '__main__':

    name = '1zoq_1'
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    cut = 2.0

    # original data
    TOP = '../../../top.prmtop'
    TRAJ = '/pod-v1/talant/Abeta/1z0q14SB/06.1.md/1z0q_1_freq10.dcd'
    traj = pt.iterload(TRAJ, TOP, stride=1)

    # DL generated synthetic data
    TRAJ_s = '../step3/' + 'SyntheticSR_Removed.dcd'
    traj_s = pt.iterload(TRAJ_s, TOP, stride=1)
    n_frames2 = traj_s.shape[0]
    if end > n_frames2:
        end = n_frames2

    
    indx_match, indx_no_match, p, count = FindCloseConfs(
        start, end, traj, traj_s, cut=cut, mask=['!@H*'])
    
    # output trajs
    pt.write_traj('out/' + name + str(start) + '_' + str(end) +
                  'SyntheticNovel.dcd', traj_s[indx_no_match], overwrite=True)
    
    print('Fraction of conformations under heavy atom cutoff' +
          str(cut) + 'similar to original MD confs', p)

    print('Number of conformations under heavy atom cutoff ' +
          str(cut) + 'similar to original MD confs', count)

