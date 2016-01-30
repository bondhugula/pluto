#!/bin/bash

export BASE=`pwd`

# ------------------------------------------------
# Utility functions
# ------------------------------------------------

# Tries to apply a patch from the given directory. Exits the program if it
# fails
check_and_apply_patch() {
  local path="$1"
  local patch="$2"

  if [ \( -d ${BASE}/$path \) -a \( ! -f ${BASE}/$path/.`basename $patch` \)  ] ; then
    echo -e "\nTrying to apply patch from directory: `pwd`"
    `git apply --check --directory=$path $patch`
    exit_status=$?

    if [ $exit_status -ne 0 ]; then
      echo -e "Failed to apply patch $patch in directory ${BASE}/$path"
      exit 1
    else
      cd ${BASE}/$path; git am $patch
      # it could have still failed (note that this is git-am, not git apply)
      exit_status=$?
      cd -

      if [ $exit_status -ne 0 ]; then
          echo -e "git-am failed on patch $patch in directory ${BASE}/$path"
          exit 1
      else
          touch ${BASE}/$path/.`basename $patch`
          echo -e "Successfully applied patch $patch in directory ${BASE}/$path"
      fi
    fi
  fi
}

# ------------------------------------------------
# Apply patches if any
# ------------------------------------------------
# check_and_apply_patch "pet" "${BASE}/patches/pet_configure_ac.patch"
check_and_apply_patch "pet" "${BASE}/patches/0001-pet-stmt-collect-accesses.patch"
check_and_apply_patch "isl" "${BASE}/patches/0001-isl-ast-build-overwrite.patch"
check_and_apply_patch "isl" "${BASE}/patches/0002-isl-dim-wise-single-valued-and-isl_pw_aff_map-functi.patch"
check_and_apply_patch "cloog-isl" "${BASE}/patches/0001-clast-multi-level-mpi-parallel-support.patch"
check_and_apply_patch "cloog-isl" "${BASE}/patches/0002-clast-data-dist-decls-print-support.patch"
check_and_apply_patch "cloog-isl" "${BASE}/patches/0003-cloog-invariant-decls-print-option.patch"
check_and_apply_patch "cloog-isl" "${BASE}/patches/0004-allow-appending-a-suffix-to-temporary-variables-for-.patch"
check_and_apply_patch "cloog-isl" "${BASE}/patches/0005-clast-introduce-lbp-ubp-as-local-variables.patch"
