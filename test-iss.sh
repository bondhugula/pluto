
#!/bin/bash

cd test
make clean > /dev/null 2>&1
cd ..

TESTS="test/jacobi-1d-periodic.c \
test/jacobi-2d-periodic.c \
test/jacobi-2d-mod.c"

for file in $TESTS; do
	echo -e "$file" 
	./polycc --pet --iss $file $*  -o test_temp_out.pluto.c
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

