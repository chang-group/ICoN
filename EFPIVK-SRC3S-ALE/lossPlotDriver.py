import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import pandas as pd

bg_cl = '#1b212c'  # Background color
tcl = 'yellow'  # Text color
legend_cl = 'white'  # Legend Color

fig = plt.figure(figsize=(10, 7), facecolor=bg_cl)
plt.rcParams['axes.facecolor'] = bg_cl
plt.tick_params(colors=tcl)
white_color = 'white'

filename = "./out/log_train_val_Cyclic_Peptide_Driver.txt"
# Check if the file exists
# if os.path.exists(filename):
#     dat = pd.read_csv(filename, sep="\t", header=None)
#     dat.columns = ["Epoch", f"LossT_", f"LossV_"]
#
#     t = dat["Epoch"].values
#     lossV = dat[f"LossV_"].values
#
#     lossT = dat[f"LossT_"].values
#
#     # Plot validation loss for each file
#     plt.semilogy(t, lossV, ls='--', lw=1.2, label=f"ValidationTorsions Error")
#     plt.semilogy(t, lossT, ls='--', lw=1.2, label=f"TrainingTorsions Error")
if os.path.exists(filename):
    dat = pd.read_csv(filename, sep="\t", header=None)
    dat.columns = ["Epoch", f"Loss_DihedralT_", f"Loss_DihedralV_", f"Loss_XYZT_", f"Loss_XYZV_"]

    t = dat["Epoch"].values
    loss1T = dat[f"Loss_DihedralT_"].values
    loss1V = dat[f"Loss_DihedralV_"].values

    loss2T = dat[f"Loss_XYZT_"].values
    loss2V = dat[f"Loss_XYZV_"].values

    # Plot validation loss for each file
    plt.semilogy(t, loss1T, ls='--', lw=1.2, label=f"Dihedral Training Loss")
    plt.semilogy(t, loss1V, ls='--', lw=1.2, label=f"Dihedral Validation Loss")
    plt.semilogy(t, loss2T, ls='--', lw=1.2, label=f"XYZ Training Loss")
    plt.semilogy(t, loss2V, ls='--', lw=1.2, label=f"XYZ Validation Loss")

plt.xlabel('Iteration', labelpad=2, color=tcl)
plt.ylabel('Loss', labelpad=2, color=tcl)
plt.legend()
plt.show()
