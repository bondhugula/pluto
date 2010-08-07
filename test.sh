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
test/tricky2.c \
test/tricky3.c \
test/tce-4index-transform.c"

for file in $TESTS; do
	echo -e "$file" 
	./polycc $file $* 
    if [ $? -ne 0 ]; then
        echo Failed!
    #else
		# Silence on success - viva la Unix
        # echo -e "Successful!"
    fi
done

echo
