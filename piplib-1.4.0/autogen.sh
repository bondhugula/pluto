#!/bin/sh
test -d autoconf || mkdir autoconf
libtoolize -c --force
aclocal -I m4
automake -a -c --foreign
autoconf
