## Instruction
1.	Go to the /output/pdb/step1 folder and make folders (i.e., rst7 and output), and within those folders, make another two folders (i.e., mdout and ncrst). 
2.	Please provide the path of the topology .prmtop and NonLinIntAll.dcd files in the GetStartingPoints.cppt. Now you can run the GetStartingPoints.cppt (It requires pytraj).
 
3.	Provide the path of the topology .prmtop file in the run_min.sh file. To run the minimization with Sander (igb=8) under AmberTools, you need to know the number of frames. Run the run_min.sh:
	./run_min.sh 1 10000 20
	First arg – initial conformation id, second is the last conformation id, and third arg – number of conformations to run in parallel. 
4.	Extract all the energy terms after minimization from .mdout files by going to /output directory and running the GetEnergires.py script.
5.	Converts all .ncrst files to a single .dcd file by using the ncrst2dcd.py script. The combined trajectories should be generated in the same directory. (Before using the script, make sure you define the path to the topology .prmtop file and ncrst folder accordingly)
Please refer to the README files in each folder for steps 2, 3, and 4.


Here we performe 

1. Split all ICoN generated conformations into individual files in `.rst7` format and save in rst7 folder
   - Inside 'GenStartingPoints.cppt' file, provide path to '.dcd' and `.prmtop`
   - `cpptraj -i GenStartingPoints.cppt`
   - It requires `pytraj`

2. Run minimization with Sander (igb=8) under AmberTools
   - `./run_min.sh 1 10000 20`
     - first arg - initial conformation id
     - second arg is last conformation id
     - 3rd arg - number of conformations to run in parallell
     
3. Extract all energy terms after minimization from '.mdout' files
   - Navigate to 'output' directory
     - `cd output/`
   - Run python code to extract energies  
      - `python GetEnergies.py`
   - make sure to change output file name, or it will use existing name and outpur `.csv` file.
4. Converts all '.ncrst' minimized conformations to single .dcd traj file
     - `python ncrst2dcd.py` 
     - All '.ncrst' files must reside inside 'ncrst' folder
     - Combined traj should be written to current dir


