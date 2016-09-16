noinst_LTLIBRARIES = libgslspmatrix.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_spmatrix.h

libgslspmatrix_la_SOURCES = spcompress.c spcopy.c spgetset.c spio.c spmatrix.c spoper.c spprop.c spswap.c

AM_CPPFLAGS = -I$(top_srcdir)

noinst_HEADERS = avl.c

TESTS = $(check_PROGRAMS)

test_LDADD = libgslspmatrix.la ../spblas/libgslspblas.la ../test/libgsltest.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../sys/libgslsys.la ../err/libgslerr.la ../utils/libutils.la ../rng/libgslrng.la

test_SOURCES = test.c
