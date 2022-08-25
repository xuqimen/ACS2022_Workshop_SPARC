#!/bin/bash

# use seq -f "%.2f" <start> <step> <end> to generate floating point sequence 
for h in `seq -f "%.2f" 0.20 .05 0.80`
do
    echo $h
    cp -r template/ mesh_$h
    sed -i "s|MESH_SPACING: 0.10|MESH_SPACING: $h|g" mesh_$h/sprc-calc.inpt 
done



