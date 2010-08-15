#! /bin/sh

libtoolize -c --force
aclocal -I m4
autoheader
automake --gnu --add-missing
autoconf

cd candl-0.4.0/
./autogen.sh
cd ..

cd clan-0.6.0/
./autogen.sh
cd ..

cd polylib-5.22.5
./autogen.sh
cd ..

cd cloog-isl
./autogen.sh
cd ..

cd piplib-1.4.0
./autogen.sh
cd ..
