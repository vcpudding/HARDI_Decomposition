#!/bin/bash
for bval in 2000
do
    for snr in 10 20 30
    do
        for denoise in 1 0
        do
            for z in {1...5}
            do
                for x in {1...10}
                do
                    for y in {1...10}
                    do
                        COMMAND=`printf "procVoxel('Training_3D_SF',%d,%d,%d,%d,%d,%d);quit" ${x} ${y} ${z} ${bval} ${snr} ${denoise}`
                        /usr/sci/local/matlab/bin/matlab -nodesktop -nosplash - nojvm -r ${COMMAND} &
                    done
                done
                wait
            done
            wait
            /usr/sci/local/matlab/bin/matlab -nodesktop -nodisplay -nosplash - nojvm -r `printf "parResults('Training_3D_SF',%d,%d,%d);quit" ${bval} ${snr} ${denoise}`
        done
    done
done
