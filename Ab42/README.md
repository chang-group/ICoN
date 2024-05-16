
<div align="center">
  <a href="https://www.youtube.com/watch?v=Yo3NVWAFruA"><img src="https://img.youtube.com/vi/Yo3NVWAFruA/0.jpg" alt
="IMAGE ALT TEXT"></a>
</div>

## ICoN - generative model for sampling highly flexible protein conformations. 

Total number of [model parameters](https://drive.google.com/file/d/17UWB6yphaCizXIb_FPon4H4b_-NThL4V/view?usp=sharing) -- 32721156

### Generation of synthetic conformations for Ab42 with ICoN model

## Contents
- Training
  - Training and validation sets for [MDrun1](https://drive.google.com/file/d/1hTl-_AWGQG7ZYrxbCy40-GyiDe6yX5uh/view?usp=sharing),
  [MDrun2](https://drive.google.com/file/d/1pkxOiUsPgQJpIh4vhQVD6sXj1nmGXPZ6/view?usp=sharing),
  [MDrun3](https://drive.google.com/file/d/1PBHn17eIpDwhHz0Go5Lyfqfx8DaWmdG8/view?usp=sharing),
  [MDrun4](https://drive.google.com/file/d/1UdiL5UhbtoS2pNQ6iuuZhdcOiKVn_v_1/view?usp=sharing).
  - Promtop file for [Ab42](https://drive.google.com/file/d/1JEtWP2Qj9CbuidBTE3MZdI7GGxJEBO4Z/view?usp=drive_link)
  - Provide path to all input files(.dcd, .prmtop) in `train.py`, then run:
  - `python train.py`
  
- Reconstruction
  - Provide path to validation set, then run:
  -`python reconstruct.py`
  
- Sampling
  - To avoid running out of GPU memory, we generate conformations in chunks (batches) by deviding entire set of synthetic conformations into 10 chunks. We provide chunk it and run it multiple times:
  - `for i in {0..9}; do python sample.py $i; done` 


