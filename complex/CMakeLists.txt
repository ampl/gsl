add_library(gslcomplex math.c inline.c)

add_gsl_test(gslcomplex-test
  test.c results.h results1.h results2.h results_real.h results_zreal.h
  LIBS gslcomplex gslerr gslieeeutils gslerr gsltest gslsys utils)
