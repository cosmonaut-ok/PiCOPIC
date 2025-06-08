#!/bin/sh -e

# generate the configure script

# touch -a install-sh

libtoolize -c -i --force
aclocal --install
autoheader -f
automake -f -a -c
autoheader -f
autoreconf -f -i
