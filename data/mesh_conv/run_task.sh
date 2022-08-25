#!/bin/bash

# use seq -f "%.2f" <start> <step> <end> to generate floating point sequence 
EXE="~/SPARC/lib/sparc"

#for i in np_*/; 
for i in mesh_*/; 
do
    cd $i;
    echo "mesh = $i"
    # qsub run.pbs
    mpirun -np 12 ${EXE} -name sprc-calc > log
    cd ..;
done

