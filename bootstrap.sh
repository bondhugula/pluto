#! /bin/sh

autoreconf -vfi

cd candl-0.4.0/
./autogen.sh
cd ..

cd clan-0.6.0/
./autogen.sh
cd scoplib/scoplib-0.2.0
autoreconf -vfi
cd ../..
cd ..

(cd isl; ./autogen.sh)

cd cloog-isl
./autogen.sh
cd ..

cd piplib
./autogen.sh
cd ..
