
<div align="center">
  <a href="https://www.youtube.com/watch?v=xchUG0sSLqI"><img src="https://img.youtube.com/vi/xchUG0sSLqI/0.jpg" alt="IMAGE ALT TEXT"></a>
</div>

## ICoN - generative model for sampling highly flexible protein conformations. 


Total number of [model parameters](https://drive.google.com/file/d/1TuqUo0TqlmM1IThc9_B4M_uDjGPHDL1m/view?usp=drive_link) --  3294964

### Generation of synthetic conformations for aB-cristalling57-69 with ICoN model

## Instruction
1.	In the train.py file, provide the path to all input files and change the trajectory and topology file names accordingly. (In train.py file, 68 and 69 lines)
2.	2.	Now you can start training your data with the `Python train.py` command line.
3.	The raw data for the loss function will be produced in /output/log_train_val.txt. The users can then check the loss function using the/visual/loss.py script.
4.	Provide the path to the validation set, which would be providing the path to the trajectory .dcd and topology .prmtop files in lines 64 and 65 in the reconstruct.py file. The reconstructed .dcd file would be generated in the /output/pdb folder.
5.	Save a .pdb (i.e., 0-frame.pdb) file from one of the frames in your original trajectory .dcd file and then transfer the 0-frame.pdb file to the /ICON_Model/TRAJ folder. 
6.	Provide the path of the 0-frame.pdb file and the topology .prmtop file in sample.py in lines 113 and 114.  The output file for sampling would be generated in /output/pdb/NonLinIntAll.dcd.
Please go to the output/pdb/step1 directory for further details of next step.


## Contents
- Training
  - [Training and validation sets](https://drive.google.com/file/d/1-VlshgKtz4Fs6p5dzLS3cv6Q-6AquoIL/view?usp=drive_link)
  - Promtop file for [aB-cristalling57-69](https://drive.google.com/file/d/10nbKLLoAYFxIaogi6aOY62uopgAbbAi5/view?usp=drive_link)
  - Provide path to all input files(.dcd, .prmtop) in `train.py`, then run:
  - `python train.py`
  
- Reconstruction
  - Provide path to validation set, then run:
  - `python reconstruct.py`
  
- Sampling
  - Sample synthetic conformations in large batch
  - `python sample.py`




