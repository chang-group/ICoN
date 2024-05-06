Here we  

1. Eliminate bad conformations (with clashes) using energy cutoff
- It requires `pytraj`
   - Inputs are:
     - '.dcd' output file of previous step, all synthetic conformations
     - '.csv' file that has all energies
   - Output:
     '.dcd' file, remaing conformations

