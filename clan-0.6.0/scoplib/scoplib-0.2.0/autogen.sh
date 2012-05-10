#! /bin/sh

aclocal -I autoconf
libtoolize --force --copy
autoreconf -vfi
