import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import pandas as pd

##=====================================================================
plot_fn= 'figs/Loss.png'
bg_cl = 'white' #'#1b212c'#'#1A1A1A'#'black'

fig = plt.figure(figsize=(5.,3.5), facecolor=bg_cl)
tcl = 'black'

dat = pd.read_csv("../output/log_train_val.txt", sep="\t", header=None)
dat.columns=["Epoch","LossT", "LossV"]
t = dat["Epoch"].values
lossT = dat["LossT"].values
lossV = dat["LossV"].values


plt.rcParams['axes.facecolor'] = bg_cl
plt.tick_params(colors=tcl)


plt.xlabel('Epoch', labelpad=2, color=tcl)
plt.ylabel('Loss', labelpad=2, color=tcl)
plt.plot(t, lossT, ls='--', lw=1.2, color='blue', alpha=0.5)
#plt.semilogx(t, lossT, ls='-', lw=1.2)
plt.plot(t, lossV, ls='-', lw=1.2, color='purple', alpha=0.5)
#plt.semilogy(t, lossV, ls='--', lw=1.2)
plt.legend(["Train", "Validation"], labelcolor = 'linecolor')
#plt.xlim(0.0, 30000)
plt.grid()

fig.set_tight_layout(True)
plt.show()
fig.savefig(plot_fn)
os.system("epscrop %s %s" % (plot_fn, plot_fn))
