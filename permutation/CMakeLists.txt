add_library(gslpermutation
  init.c file.c permutation.c permute.c canonical.c inline.c)

add_gsl_test(gslpermutation-test test.c
  LIBS gslpermutation gslvector gslblock gslieeeutils
  gslerr gsltest gslsys utils)
