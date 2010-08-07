#! /bin/sh

aclocal
automake --gnu --add-missing
autoconf
autoheader

cd candl-0.4.0/
./autogen.sh
cd ..

cd clan-0.6.0/
./autogen.sh
cd ..

cd polylib-5.22.3
./autogen.sh
cd ..

cd cloog-isl
./autogen.sh
cd ..

cd piplib-1.4.0
./autogen.sh
cd ..
