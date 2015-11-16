noinst_LTLIBRARIES = libgslrstat.la

pkginclude_HEADERS = gsl_rstat.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslrstat_la_SOURCES = rstat.c rquantile.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslrstat.la ../statistics/libgslstatistics.la ../sort/libgslsort.la ../ieee-utils/libgslieeeutils.la ../randist/libgslrandist.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../complex/libgslcomplex.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la ../vector/libgslvector.la


