sleep `expr 3 - $1`
echo $1
if [ -d test$1 ]; then
    rmdir test$1
fi
mkdir test$1
