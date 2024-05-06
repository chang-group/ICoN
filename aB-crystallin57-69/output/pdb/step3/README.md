
Here we eliminate self repeats based on 1A heavy atom rmsd cutoff
- It requires `pytraj`
- Inputs are:
  - '.dcd' output file of previous step & '.prmtop' file
- Output:
  - '.dcd' file, remaing conformations
- Make sure all input and output path are correct then run
  - `python RmSelfRepeats.py`


