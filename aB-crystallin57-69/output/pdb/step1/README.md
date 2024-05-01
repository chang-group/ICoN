Here we performe 

1. Split all ICoN generated conformations into individual files in `.rst7` format and save in rst7 folder
   - Inside 'GenStartingPoints.cppt' file, provide path to '.dcd' and `.prmtop`
   - `cpptraj -i GenStartingPoints.cppt`

2. Run minimization with Sander (igb=8) under AmberTools
   - `./run_min.sh 1 10000 20`
     - first arg - initial conformation id
     - second arg is last conformation id
     - 3rd arg - number of conformations to run in parallell
     
3. Extract all energy terms after minimization from '.mdout' files
   - `python GetEnergies.py`
     - make sure to change output file name, or it will use existing name and outpur `.csv` file.
4. Converts all '.ncrst' minimized conformations to single .dcd traj file
     - `python ncrst2dcd.py` 



