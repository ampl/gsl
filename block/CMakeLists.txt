add_library(gslblock init.c file.c block.c)

add_gsl_test(gslblock-test test.c
  LIBS gslblock gslieeeutils gslerr gsltest gslsys utils)
