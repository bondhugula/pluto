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
    echo "Failed test case $file!"
    exit_value=${test_exit_value}
  fi
}


#
# MAIN
#

# Tests to run
TESTS="test/2d-bidirec.c \
  test/reverse-iss.c"

# Add trap for clean
trap cleanup SIGINT exit

# Make tests
cd test || exit 1
make clean > /dev/null 2>&1
cd .. || exit 1

# Run tests
exit_value=0

OPTS="--iss --islsolve"
for file in $TESTS; do
  cmd="./polycc $OPTS $file $* -o test_temp_out.pluto.c"
  run_test "$file" "$cmd"
done

OPTS="--iss --islsolve --lbtile --parallel"
for file in $TESTS; do
  cmd="./polycc $OPTS $file $* -o test_temp_out.pluto.c"
  run_test "$file" "$cmd"
done

# End
echo 
if [ "${exit_value}" -ne 0 ]; then
  echo "ERROR: Some test failed. Please check errors above"
fi
exit ${exit_value}

