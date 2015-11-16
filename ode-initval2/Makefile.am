noinst_LTLIBRARIES = libgslodeiv2.la 

pkginclude_HEADERS = gsl_odeiv2.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslodeiv2_la_SOURCES = control.c cstd.c cscal.c evolve.c step.c rk2.c rk2imp.c rk4.c rk4imp.c rkf45.c rk8pd.c rkck.c bsimp.c rk1imp.c msadams.c msbdf.c driver.c

noinst_HEADERS = odeiv_util.h step_utils.c rksubs.c modnewton1.c control_utils.c

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_LDADD = libgslodeiv2.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../permutation/libgslpermutation.la ../vector/libgslvector.la ../block/libgslblock.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la 

test_SOURCES = test.c

