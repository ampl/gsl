noinst_LTLIBRARIES = libgslmultilarge.la 

pkginclude_HEADERS = gsl_multilarge.h

libgslmultilarge_la_SOURCES = multilarge.c normal.c tsqr.c

AM_CPPFLAGS = -I$(top_srcdir)

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslmultilarge.la ../test/libgsltest.la ../multifit/libgslmultifit.la ../eigen/libgsleigen.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../permutation/libgslpermutation.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../sys/libgslsys.la ../utils/libutils.la ../rng/libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../complex/libgslcomplex.la ../min/libgslmin.la
