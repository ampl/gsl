noinst_LTLIBRARIES = libgslmultifit_nlinear.la 

pkginclude_HEADERS = gsl_multifit_nlinear.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslmultifit_nlinear_la_SOURCES = cholesky.c convergence.c covar.c dogleg.c fdf.c fdfvv.c fdjac.c lm.c mcholesky.c qr.c scaling.c subspace2D.c svd.c trust.c

noinst_HEADERS =        \
common.c                \
nielsen.c               \
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
test_exp1.c             \
test_fdf.c              \
test_gaussian.c         \
test_hahn1.c            \
test_helical.c          \
test_jennrich.c         \
test_kirby2.c           \
test_kowalik.c          \
test_lin1.c             \
test_lin2.c             \
test_lin3.c             \
test_meyer.c            \
test_meyerscal.c        \
test_osborne.c          \
test_penalty1.c         \
test_penalty2.c         \
test_powell1.c          \
test_powell2.c          \
test_powell3.c          \
test_rat42.c            \
test_rat43.c            \
test_rosenbrock.c       \
test_rosenbrocke.c      \
test_roth.c             \
test_thurber.c          \
test_vardim.c           \
test_watson.c           \
test_wnlin.c            \
test_wood.c

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslmultifit_nlinear.la ../eigen/libgsleigen.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../sort/libgslsort.la ../statistics/libgslstatistics.la ../vector/libgslvector.la ../block/libgslblock.la  ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../poly/libgslpoly.la
