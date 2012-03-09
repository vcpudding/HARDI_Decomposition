BVAL=$1
SNR=$2
DENOISE=$3
ZS=($(seq 1 5))
ROWS=($(seq 1 16))
COLS=($(seq 1 16))
COMMAND_FORMAT="procVoxel('Training_3D_SF',%d,%d,%d,%d,%d,%d);"
for z in ${ZS[@]}; do
    wait
    for x in ${COLS[@]}; do
        for y in ${ROWS[@]}; do
            COMMAND=`printf ${COMMAND_FORMAT} ${x} ${y} ${z} ${BVAL} ${SNR} ${DENOISE}`
            echo ${COMMAND}
            /usr/sci/local/matlab/bin/matlab -nodesktop -nosplash - nojvm -r ${COMMAND} &
        done
    done
done
