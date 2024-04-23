



## ICoN - generative model for sampling highly flexible protein conformations. 

<video src='visual/ab13.mp4'  width=200/>

Total number of [model parameters](https://drive.google.com/file/d/1TuqUo0TqlmM1IThc9_B4M_uDjGPHDL1m/view?usp=drive_link) --  3294964

### Generation of synthetic conformations for aB-cristalling57-69 with ICoN model

## Contents
- Training
  - Provide path to all input files(.dcd, .prmtop) in `train.py`, then run:
  - `python train.py`
  
- Reconstruction
  - Provide path to validation set, then run:
  - `python reconstruct.py`
  
- Sampling
  - Sample synthetic conformations in large batch
  - `python sample.py` 


