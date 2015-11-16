#! /bin/sh

# Run this to generate all the auto-generated files needed by the GNU
# configure program
#libtoolize --automake
#aclocal
#autoheader
#automake --add-missing --gnu --force-missing
#autoconf
autoreconf -i -f -v
echo "Now use ./configure --enable-maintainer-mode"
