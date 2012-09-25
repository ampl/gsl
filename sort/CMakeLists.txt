add_library(gslsort
  sort.c sortind.c sortvec.c sortvecind.c subset.c subsetind.c)

add_gsl_test(gslsort-test test.c
  LIBS gslsort gslpermutation gslvector gslblock
  gslieeeutils gslerr gsltest gslsys)
