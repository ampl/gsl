noinst_LTLIBRARIES = libgslbst.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_bst.h gsl_bst_types.h gsl_bst_avl.h gsl_bst_rb.h

libgslbst_la_SOURCES = avl.c bst.c rb.c trav.c

AM_CPPFLAGS = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslbst.la ../test/libgsltest.la ../sys/libgslsys.la ../err/libgslerr.la ../sort/libgslsort.la ../rng/libgslrng.la
