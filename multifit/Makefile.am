noinst_LTLIBRARIES = libgslmultifit.la 

pkginclude_HEADERS = gsl_multifit.h gsl_multifit_nlin.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslmultifit_la_SOURCES = gcv.c multilinear.c multiwlinear.c work.c lmniel.c lmder.c fsolver.c fdfsolver.c fdfridge.c fdjac.c convergence.c gradient.c covar.c multirobust.c robust_wfun.c multireg.c

noinst_HEADERS =        \
linear_common.c         \
lmutil.c                \
lmpar.c                 \
lmset.c                 \
lmiterate.c             \
lmmisc.c                \
qrsolv.c                \
test_bard.c             \
test_beale.c            \
test_biggs.c            \
test_box.c              \
test_boxbod.c           \
test_brown1.c           \
test_brown2.c           \
test_brown3.c           \
test_eckerle.c          \
test_enso.c             \
test_estimator.c        \
test_exp1.c             \
test_filip.c            \
test_gaussian.c         \
test_hahn1.c            \
test_helical.c          \
test_jennrich.c         \
test_kirby2.c           \
test_kowalik.c          \
test_lin1.c             \
test_lin2.c             \
test_lin3.c             \
test_linear.c           \
test_longley.c          \
test_meyer.c            \
test_meyerscal.c        \
test_nelson.c           \
test_nonlinear.c        \
test_osborne.c          \
test_penalty1.c         \
test_penalty2.c         \
test_pontius.c          \
test_powell1.c          \
test_powell2.c          \
test_powell3.c          \
test_rat42.c            \
test_rat43.c            \
test_reg.c              \
test_rosenbrock.c       \
test_rosenbrocke.c      \
test_roth.c             \
test_shaw.c             \
test_thurber.c          \
test_vardim.c           \
test_watson.c           \
test_wnlin.c            \
test_wood.c

check_PROGRAMS = test #demo

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslmultifit.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../sort/libgslsort.la ../statistics/libgslstatistics.la ../vector/libgslvector.la ../block/libgslblock.la  ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../min/libgslmin.la

#demo_SOURCES = demo.c
#demo_LDADD = libgslmultifit.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../randist/libgslrandist.la ../rng/libgslrng.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la

