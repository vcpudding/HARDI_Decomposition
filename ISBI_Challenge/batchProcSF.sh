#!/bin/bash
BVAL=$1
SNR=$2
DENOISE=$3
Z=1
ROWS=($(seq 1 2))
COLS=($(seq 1 2))
COMMAND_FORMAT="procVoxel('Training_SF',%d,%d,%d,%d,%d,%d);quit;"
i=1
#for x in ${COLS[@]}; do
#    for y in ${ROWS[@]}; do
for x in {1..2}
do
    for y in {1..2}
    do
        COMMAND=`printf ${COMMAND_FORMAT} ${x} ${y} ${Z} ${BVAL} ${SNR} ${DENOISE}`
        #echo ${COMMAND}
        /usr/sci/local/matlab/bin/matlab -nodesktop -nodisplay -nosplash - nojvm -r ${COMMAND} &
        echo $i
        i=`expr $i + 1`
    done
done

