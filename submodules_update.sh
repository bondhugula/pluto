#!/bin/sh

set -e

submodule_folder_list=$(grep "path" .gitmodules | awk '{ print $NF }')

for submodule in ${submodule_folder_list}; do
  echo "*** Update submodule $submodule ***"
  (cd "$submodule"; git pull origin master)
  git add "$submodule"
  if git commit -m "Update submodule $submodule"; then
    echo "*** Commit changes on submodule $submodule ***"
  fi
done
