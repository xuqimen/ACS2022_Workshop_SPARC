#!/bin/bash
for i in np_*/; 
do
    cd $i;
    echo $i
    qsub run.pbs
    cd ..;
done

