for i in 1 2 3
do
    lastidx=`expr ${i} - 1`
    lastfile=test${lastidx}/done
    while [ ! -d test$lastidx ]
    do
        sleep 1
    done
    sh batchTest.sh ${i} &
done
