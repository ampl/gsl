noinst_LTLIBRARIES = libgslmovstat.la

pkginclude_HEADERS = gsl_movstat.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslmovstat_la_SOURCES = \
  alloc.c                  \
  apply.c                  \
  fill.c                   \
  funcacc.c                \
	madacc.c                 \
	medacc.c                 \
	mmacc.c                  \
	movmad.c                 \
	movmean.c                \
	movmedian.c              \
	movminmax.c              \
	movsum.c                 \
	movSn.c                  \
	movQn.c                  \
	movqqr.c                 \
	movvariance.c            \
	mvacc.c                  \
	qnacc.c                  \
	qqracc.c                 \
	snacc.c                  \
	sumacc.c

noinst_HEADERS = deque.c ringbuf.c test_mad.c test_mean.c test_median.c test_minmax.c test_Qn.c test_qqr.c test_Sn.c test_sum.c test_variance.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslmovstat.la ../statistics/libgslstatistics.la ../sort/libgslsort.la ../ieee-utils/libgslieeeutils.la ../randist/libgslrandist.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../complex/libgslcomplex.la ../err/libgslerr.la ../test/libgsltest.la ../vector/libgslvector.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../block/libgslblock.la ../sys/libgslsys.la ../utils/libutils.la
