
<p align="center">
  <img src="docs/figs/ICoN.png">
</p>

## ICoN model for sampling highly flexible protein conformations. 

1. ICoN is trained on ~10000 of fully atomistic and highly flexible conformations

2. It uses vector internal coordinate representation as input features- vBAT

3. It can train ~10K conformations in a few mins, and generate ~100K conformations in less than a min.




## Dependancies
- `python` 3.8 >
  - pytorch
  - `MDAnalysis`
  

## Installation

Clone the repository:

*git clone https://github.com/chang-group/BondAngleTorsion.git*

then run the following command:

*python train.py*

*python reconstruct.py*

*python sample.py*



