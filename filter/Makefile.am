noinst_LTLIBRARIES = libgslfilter.la

pkginclude_HEADERS = gsl_filter.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslfilter_la_SOURCES = gaussian.c impulse.c median.c rmedian.c

noinst_HEADERS = ringbuf.c test_impulse.c test_gaussian.c test_median.c test_rmedian.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslfilter.la ../movstat/libgslmovstat.la ../statistics/libgslstatistics.la ../sort/libgslsort.la ../ieee-utils/libgslieeeutils.la ../randist/libgslrandist.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../complex/libgslcomplex.la ../err/libgslerr.la ../test/libgsltest.la ../vector/libgslvector.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../block/libgslblock.la ../poly/libgslpoly.la ../sys/libgslsys.la ../utils/libutils.la
