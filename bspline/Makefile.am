noinst_LTLIBRARIES = libgslbspline.la 

pkginclude_HEADERS = gsl_bspline.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslbspline_la_SOURCES = bspline.c eval.c gram.c greville.c integ.c interp.c ls.c old.c inline.c

noinst_HEADERS = test_gram.c test_greville.c test_integ.c test_interp.c test_ls.c

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_LDADD = libgslbspline.la ../multifit/libgslmultifit.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../complex/libgslcomplex.la ../cblas/libgslcblas.la ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la ../statistics/libgslstatistics.la ../poly/libgslpoly.la ../integration/libgslintegration.la ../rng/libgslrng.la ../sort/libgslsort.la

test_SOURCES = test.c
