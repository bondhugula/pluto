#! /bin/sh

libtoolize --force --copy
aclocal -I autoconf
automake -a -c --foreign
autoconf -f
