add_library(gslieeeutils print.c make_rep.c gsl_ieee_utils.h env.c fp.c read.c)

add_gsl_test(gslieeeutils-test test.c
  LIBS gslieeeutils gslerr gsltest gslsys utils)
