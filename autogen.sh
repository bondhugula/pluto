#!/usr/bin/env bash

export BASE=`pwd`

# ------------------------------------------------
# (Re)Generate autotools files
# ------------------------------------------------

echo -e "\n*** Running autotools on pluto ***"
autoreconf -vi

echo -e "\n*** Running autotools on isl ***"
(cd isl; ./autogen.sh)

echo -e "\n*** Running autotools on cloog-isl ***"
cd cloog-isl
./autogen.sh
cd ..

echo -e "\n*** Running autotools on piplib ***"
cd piplib
./autogen.sh
cd ..

echo -e "\n*** Running autotools on polylib ***"
cd polylib
autoreconf -vi
cd ..

echo -e "\n*** Running autotools on openscop ***"
cd openscop
autoreconf -vi
cd ..

echo -e "\n*** Running autotools on candl ***"
cd candl
autoreconf -vi
cd ..

echo -e "\n*** Running autotools on clan ***"
cd clan
./autogen.sh
cd ..

cd ..
