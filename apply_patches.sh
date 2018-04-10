#!/bin/bash

BASE=$(pwd)

# ------------------------------------------------
# Utility functions
# ------------------------------------------------

# Tries to apply a patch from the given directory. Exits the program if it
# fails
check_and_apply_patch() {
  local path
  local patch
  local path_complete
  local patch_basename
  local path_to_patch_file

  path=$1
  patch=$2
  path_complete=${BASE}/${path}
  patch_basename=$(basename "${patch}")
  path_to_patch_file=${path_complete}/.${patch_basename}

  if [ -d "${path_complete}" ] && [ ! -f "${path_to_patch_file}" ]; then
    echo -e "\nTrying to apply patch from directory: $BASE"
    git apply --check --directory="$path" "$patch"
    exit_status=$?

    if [ ${exit_status} -ne 0 ]; then
      echo -e "Failed to apply patch ${patch} in directory ${BASE}/${path}"
      exit ${exit_status}
    else
      cd "${BASE}"/"${path}" || exit 1
      git am "$patch"
      # It could have still failed (note that this is git-am, not git apply)
      exit_status=$?
      cd - || exit 1

      if [ ${exit_status} -ne 0 ]; then
          echo -e "git-am failed on patch ${patch} in directory ${BASE}/${path}"
          exit ${exit_status}
      else
          touch "${path_to_patch_file}"
          echo -e "Successfully applied patch ${patch} in directory ${BASE}/${path}"
      fi
    fi
  fi
}

check_and_apply_patch "isl" "${BASE}/patches/0001-isl-dim-wise-single_valued-and-isl_pw_aff_map-functi.patch"
