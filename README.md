
<p align="center">
  <img src="docs/figs/ICoN.png">
</p>

## Internal Coordinate Net (ICoN) for Sampling Conformational Ensembles of Highly Dynamic Proteins via Generative Deep Learning

1. ICoN is trained on ~10000 of fully atomistic and highly flexible conformations

2. It uses vector internal coordinate representation as input features- vBAT

3. It can train ~10K conformations in a few mins, and generate ~100K conformations in less than a min.



## Dependancies
- `python` 3.8 >
  - pytorch - for Deep Learning
  - `MDAnalysis` - for trajectory I/O
  - `pytraj` - for analysis

## Installation

- Clone the repository:
  *git clone https://github.com/chang-group/ICoN.git*
- Navigate to corresponding folder for further instruction
## Instuctions for users

- If users want to use this package for their system, they can download the aB-crystallin57-69 folder containing three different folders (output, src, and visual). Please follow these steps:

1.	Make a folder (i.e., ICON_Model) and copy the downloaded folder (aB-crystallin57-69) into it. (you can change the name of both folders accordingly)
2.	Make a folder (i.e., TRAJ) within the ICON_Model folder to copy the trajectory .dcd and topology files .prmtop.
Please go to the aB-crystallin57-69 folder for further details of each step.
![image](https://github.com/user-attachments/assets/b13314f9-7337-45fc-9d2e-6cef499ed915)



