#!/usr/bin/env bash

export BASE=`pwd`

# ------------------------------------------------
# Utility functions
# ------------------------------------------------

# Tries to apply a patch from the given directory. Exits the program if it
# fails
check_and_apply_patch() {
  local path="$1"
  local patch="$2"

  if [ \( -d $path \) -a \( ! -f $path/.`basename $patch` \)  ] ; then
    echo -e "\nTrying to apply patch from directory: `pwd`"
    `git apply --check --directory=$path $patch`
    exit_status=$?

    if [ $exit_status -ne 0 ]; then
      echo -e "Failed to apply patch $patch in directory $path"
      exit 1
    else
      `git apply --directory=$path $patch`
      touch $path/.`basename $patch`
      echo -e "Successfully applied patch $patch in directory $path"
    fi
  fi
}


# ------------------------------------------------
# Apply patches if any
# ------------------------------------------------

check_and_apply_patch "${BASE}/pet" "${BASE}/patches/pet_stmt_text.patch"
check_and_apply_patch "${BASE}/isl" "${BASE}/patches/isl_ast_build_overwrite.patch"
