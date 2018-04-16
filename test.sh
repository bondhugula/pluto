#!/bin/bash

#
# Helper functions
#

cleanup() {
  rm -f test_temp_out.pluto.c
  rm -f test_temp_out.pluto.pluto.cloog
}

# Runs command cmd and updates global exit_value on error
run_test() {
  local file=$1
  local cmd=$2
  local test_exit_value

  # Log information
  echo "$file"
  echo "$cmd"

  # Execute test
  $cmd

  # Check exit value
  test_exit_value=$?
  if [ ${test_exit_value} -ne 0 ]; then                                                                                                                                                                                                                                    
    echo -e "\e[31mFailed\e[0m" " $file"!
    exit_value=${test_exit_value}
  else
    echo -e "\e[32mPassed\e[0m"
  fi
}

  
#
# MAIN
#

# TESTS to run
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

  
# Add trap for clean
trap cleanup SIGINT exit

# Make tests
cd test || exit 1
make clean > /dev/null 2>&1
cd .. || exit 1

# Run tests
exit_value=0
for file in $TESTS; do
  cmd="./polycc $file $* -o test_temp_out.pluto.c"
  run_test "$file" "$cmd"
done

# End
echo 
if [ "${exit_value}" -ne 0 ]; then
  echo "ERROR: Some test failed. Please check errors above"
fi
exit ${exit_value}

