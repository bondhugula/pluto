#!/bin/bash

TESTS=test/heat-2d.c

for file in $TESTS; do
	echo -ne "$file "
	./src/pluto $file $* --pet --notile --noparallel  -o test_temp_out.pluto.c | FileCheck $file
    if [ $? -ne 0 ]; then
        echo -e "[\e[31mFailed\e[0m]" " $file"!
    else
        echo -e "[\e[32mPassed\e[0m]"
    fi
done

cleanup()
{
rm -f test_temp_out.pluto.c
rm -f test_temp_out.pluto.pluto.cloog
}

echo

trap cleanup SIGINT exit
