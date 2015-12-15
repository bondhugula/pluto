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
      echo Hello
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

check_and_apply_patch "isl" "${BASE}/patches/0001-isl-dim-wise-single_valued-and-isl_pw_aff_map-functi.patch"
