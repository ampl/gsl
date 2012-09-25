add_library(gslvector
  init.c file.c vector.c copy.c swap.c prop.c minmax.c oper.c
  reim.c subvector.c view.c view.h)

set(LIBS gslvector gslblock gslieeeutils gslerr gsltest gslsys utils)
add_gsl_test(gslvector-test test.c LIBS ${LIBS})
add_gsl_test(gslvector-test-static test_static.c LIBS ${LIBS})
