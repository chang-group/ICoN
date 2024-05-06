
Here we eliminate comformations by comparing to raw MD using 2A heavy atom RMSD cutoff

- It requires `pytraj`
- Inputs are:
  - '.dcd' output file of previous step & '.prmtop' file
  - also provide path to raw MD '.dcd' file
- Output:
  - '.dcd' file, remaing conformations
- Make sure all input and output path are correct then run
  - `python GetNovelConfs.py start end`
    - In order to speedup calculation one can run separate chunks in separte CPU with interval
    - start - is starting synthetic conformations ID
    - last - last synthetic conformations ID of the chunk
- Once selected conformations are outputted in chunks, one can merge them using
  - `python MergeDCD.py`