#!/usr/bin/env bash

export BASE=`pwd`

# ------------------------------------------------
# Apply patches if any
# ------------------------------------------------

cd pet
patch -p1 < "${BASE}/patches/pet_stmt_text.patch"
cd ..

cd isl/
patch -p1 < "${BASE}/patches/isl_ast_build_overwrite.patch"
cd ..
