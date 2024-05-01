#!/bin/bash
START=$1                #starting traj id, 1nd cmd arg
END=$2                  #ending traj id, 2rd cmd arg 
BATCH=$3                #number of confs to minimize in different cpus (parallell)
#GPU_ID=4 #$4           #which gpu to use
#export CUDA_VISIBLE_DEVICES=$GPU_ID

for (( i=$START; i<=$END; i++ ))
do
    sander -O -i MD.mdin \
	 -p ../top.prmtop \
	 -o output/mdout/min_fr_$i.mdout \
	 -c rst7/MD.rst7.$i \
	 -r output/ncrst/Min_$i.ncrst &  # & makes jobs concurrent, and runs in bg
                                         # make sure number of jobs less than number of cpu cores
    # runs n=BATCH jobs simulataneously using different cpus
    if (( $i % $BATCH == 0))
    then
	echo 'Waiting for batch' $i $[i+$BATCH]
	wait
    fi		          
done
