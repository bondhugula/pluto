#!/bin/sh
libtoolize -c
aclocal -I m4
autoheader
automake -a -c --foreign
autoconf
if test -f piplib/autogen.sh; then
	(cd piplib; ./autogen.sh)
fi
