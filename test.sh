#!/bin/bash

cd test
make clean > /dev/null 2>&1
cd ..

TESTS="test/jacobi-1d-imper.c \
test/jacobi-2d-imper.c \
test/matmul.c \
test/costfunc.c \
test/fdtd-2d.c \
test/heat-3d-imperfect.c \
test/seq.c \
test/gemver.c \
test/seidel.c \
test/mvt.c \
test/mxv.c \
test/mxv-seq.c \
test/mxv-seq3.c \
test/matmul-seq.c \
test/matmul-seq3.c \
test/darte.c \
test/doitgen.c \
test/polynomial.c \
test/1dloop-invar.c \
test/nodep.c \
test/simple.c \
test/fusion1.c \
test/fusion2.c \
test/fusion3.c \
test/fusion4.c \
test/fusion5.c \
test/fusion6.c \
test/fusion7.c \
test/fusion8.c \
test/fusion9.c \
test/fusion10.c \
test/negparam.c \
test/tricky1.c \
test/tricky2.c \
test/tricky3.c \
test/tricky4.c \
test/multi-stmt-lazy-lin-ind.c \
test/ludcmp.c \
test/tce-4index-transform.c \
test/noloop.c"

for file in $TESTS; do
	echo -e "$file"
	./polycc $file $*  -o test_temp_out.pluto.c
    if [ $? -ne 0 ]; then
        echo -e "\e[31mFailed\e[0m" " $file"!
    else
        echo -e "\e[32mPassed\e[0m"
    fi
done

cleanup()
{
rm -f test_temp_out.pluto.c
rm -f test_temp_out.pluto.pluto.cloog
}

echo

trap cleanup SIGINT exit
