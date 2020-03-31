noinst_LTLIBRARIES = libgslspblas.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_spblas.h

libgslspblas_la_SOURCES = spdgemm.c spdgemv.c

AM_CPPFLAGS = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_LDADD = libgslspblas.la ../spmatrix/libgslspmatrix.la ../bst/libgslbst.la ../test/libgsltest.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../sys/libgslsys.la ../err/libgslerr.la ../utils/libutils.la ../rng/libgslrng.la

test_SOURCES = test.c
