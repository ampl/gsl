noinst_LTLIBRARIES = libgslsplinalg.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_splinalg.h

libgslsplinalg_la_SOURCES = itersolve.c gmres.c

AM_CPPFLAGS = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_LDADD = libgslsplinalg.la ../spmatrix/libgslspmatrix.la ../spblas/libgslspblas.la ../bst/libgslbst.la ../test/libgsltest.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../sys/libgslsys.la ../utils/libutils.la ../rng/libgslrng.la ../err/libgslerr.la

test_SOURCES = test.c
