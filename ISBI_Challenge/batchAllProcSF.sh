#!/bin/bash
COMMAND_FORMAT="procVoxel('Training_SF',%d,%d,1,%d,%d,%d);quit;"
for bval in 3000 2000
do
    for snr in 10 20 30
    do
        for denoise in 1 #0
        do
            for x in {1..10}
            do
                for y in {1..10}
                do
                    COMMAND=`printf ${COMMAND_FORMAT} ${x} ${y} ${bval} ${snr} ${denoise}`
                    #echo ${COMMAND}
                    /usr/sci/local/matlab/bin/matlab -nodesktop -nodisplay -nosplash - nojvm -r ${COMMAND} &
                done
            done
            wait
            /usr/sci/local/matlab/bin/matlab -nodesktop -nodisplay -nosplash - nojvm -r `printf "parResults('Training_SF',%d,%d,%d);quit" ${bval} ${snr} ${denoise}`
        done
        #wait
    done
    #wait
done
