#! /bin/sh

autoreconf -vi

(cd isl; ./autogen.sh)

cd cloog-isl
./autogen.sh
cd ..

cd piplib
./autogen.sh
cd ..

cd polylib
autoreconf -vi
cd ..

cd openscop
autoreconf -vi
cd ..

cd candl
autoreconf -vi
cd ..

cd clan
#./get_submodules
#cd osl
./autogen.sh
cd ..

autoreconf -vi
cd ..
