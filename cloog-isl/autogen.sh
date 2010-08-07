#!/bin/sh
libtoolize -c --force
aclocal -I m4
autoheader
automake -a -c --foreign
autoconf
if test -f isl/autogen.sh; then
	(cd isl; ./autogen.sh)
fi
