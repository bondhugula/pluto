#!/bin/bash

cd test
make clean > /dev/null 2>&1
cd ..

TESTS="test/jacobi-1d-imper.c \
test/jacobi-2d-imper.c \
test/costfunc.c \
test/fdtd-2d.c \
test/gemver.c \
test/matmul.c \
test/seidel.c \
test/mvt.c \
test/darte.c \
test/doitgen.c \
test/polynomial.c \
test/1dloop-invar.c \
test/nodep.c \
test/simple.c \
test/test7.c \
test/test8.c \
test/tricky1.c \
test/tricky2.c \
test/tricky3.c \
test/tce-4index-transform.c"

for file in $TESTS; do
	echo -e "$file" 
	./polycc $file $*  -o test_temp_out.pluto.c
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
