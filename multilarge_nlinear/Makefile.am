noinst_LTLIBRARIES = libgslmultilarge_nlinear.la 

pkginclude_HEADERS = gsl_multilarge_nlinear.h

libgslmultilarge_nlinear_la_SOURCES = cgst.c cholesky.c convergence.c dogleg.c dummy.c fdf.c lm.c mcholesky.c scaling.c subspace2D.c trust.c

AM_CPPFLAGS = -I$(top_srcdir)

noinst_HEADERS =     \
common.c             \
nielsen.c            \
test_fdf.c           \
test_bard.c          \
test_beale.c         \
test_biggs.c         \
test_box.c           \
test_boxbod.c        \
test_brown1.c        \
test_brown2.c        \
test_brown3.c        \
test_eckerle.c       \
test_enso.c          \
test_exp1.c          \
test_hahn1.c         \
test_helical.c       \
test_gaussian.c      \
test_jennrich.c      \
test_kirby2.c        \
test_kowalik.c       \
test_lin1.c          \
test_lin2.c          \
test_lin3.c          \
test_meyer.c         \
test_meyerscal.c     \
test_osborne.c       \
test_penalty1.c      \
test_penalty2.c      \
test_powell1.c       \
test_powell2.c       \
test_powell3.c       \
test_rat42.c         \
test_rat43.c         \
test_rosenbrock.c    \
test_rosenbrocke.c   \
test_roth.c          \
test_thurber.c       \
test_vardim.c        \
test_watson.c        \
test_wood.c          \
test_wnlin.c

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslmultilarge_nlinear.la ../test/libgsltest.la ../multilarge/libgslmultilarge.la ../multifit/libgslmultifit.la ../eigen/libgsleigen.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../sys/libgslsys.la ../utils/libutils.la ../rng/libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../complex/libgslcomplex.la ../poly/libgslpoly.la
