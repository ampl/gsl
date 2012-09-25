add_library(gslmatrix
  init.c matrix.c file.c rowcol.c swap.c copy.c minmax.c prop.c
  oper.c getset.c view.c submatrix.c view.h)

set(LIBS gslmatrix gslvector gslblock gslieeeutils gslerr gsltest gslsys utils)
add_gsl_test(gslmatrix-test test.c LIBS ${LIBS})
add_gsl_test(gslmatrix-test-static test_static.c LIBS ${LIBS})
