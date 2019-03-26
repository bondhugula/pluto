
#!/bin/bash

cd test
make clean > /dev/null 2>&1
cd ..

TESTS="test/2d-bidirec.c \
test/reverse-iss.c
"

OPTS="--iss --islsolve"

for file in $TESTS; do
	echo -e "$file" 
    echo ./polycc $OPTS $file $*  -o test_temp_out.pluto.c
    ./polycc $OPTS $file $*  -o test_temp_out.pluto.c
    if [ $? -ne 0 ]; then
        echo Failed test case "$file"!
        break
    fi
done

OPTS="--iss --islsolve --parallel"
for file in $TESTS; do
	echo -e "$file" 
    echo ./polycc $OPTS $file $*  -o test_temp_out.pluto.c
    ./polycc $OPTS $file $*  -o test_temp_out.pluto.c
    if [ $? -ne 0 ]; then
        echo Failed test case "$file"!
        break
    fi
done


cleanup()
{
rm -f test_temp_out.pluto.c
rm -f test_temp_out.pluto.pluto.cloog
}

echo

trap cleanup SIGINT exit

