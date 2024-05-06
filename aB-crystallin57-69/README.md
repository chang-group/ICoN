
<div align="center">
  <a href="https://www.youtube.com/watch?v=xchUG0sSLqI"><img src="https://img.youtube.com/vi/xchUG0sSLqI/0.jpg" alt="IMAGE ALT TEXT"></a>
</div>

## ICoN - generative model for sampling highly flexible protein conformations. 


Total number of [model parameters](https://drive.google.com/file/d/1TuqUo0TqlmM1IThc9_B4M_uDjGPHDL1m/view?usp=drive_link) --  3294964

### Generation of synthetic conformations for aB-cristalling57-69 with ICoN model

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


