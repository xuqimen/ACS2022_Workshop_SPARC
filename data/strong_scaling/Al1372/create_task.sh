#!/bin/bash

# use seq -f "%.2f" <start> <step> <end> to generate floating point sequence 
#for i in 11 12 14 15 18 22 27 36 54
#for h in `seq -f "%.2f" 0.1 .05 0.9`

#for i in `seq -f "%.0f" 50 10 300`
for i in 24 48 96 192 312
do
    echo $i
    cp -r template/ np_$i
    nnode=`python3 -c "import math; print(int(math.ceil(${i}/24.0)))"` # find number of nodes
    sed -i "s|__nnode__|${nnode}|g" np_$i/run.pbs
    sed -i "s|nproc=24|nproc=$i|g" np_$i/run.pbs
    # sed -i "s|NSTATES: 1|NSTATES: $i|g" np_$i/sprc-calc.inpt 
done



