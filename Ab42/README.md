
## ICoN - generative model for sampling highly flexible protein conformations. 

Total number of [model parameters](https://drive.google.com/file/d/17UWB6yphaCizXIb_FPon4H4b_-NThL4V/view?usp=sharing) -- 32721156

### Generation of synthetic conformations for Ab42 with ICoN model

## Contents
- Training
  - Provide path to all input files(.dcd, .prmtop) in `train.py`, then run:
  - `python train.py`
  
- Reconstruction
  - Provide path to validation set, then run:
  -`python reconstruct.py`
  
- Sampling
  - To avoid running out of GPU memory, we generate conformations in chunks (batches) by deviding entire set of synthetic conformations into 10 chunks. We provide chunk it and run it multiple times:
  - `for i in {0..9}; do python sample.py $i; done` 


