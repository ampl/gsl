add_library(gsleigen
  jacobi.c symm.c symmv.c nonsymm.c nonsymmv.c herm.c hermv.c gensymm.c
  gensymmv.c genherm.c genhermv.c gen.c genv.c sort.c francis.c schur.c)

add_gsl_test(gsleigen-test test.c
  LIBS gsleigen gsltest gsllinalg gslpermutation gslblas gslcblas gslmatrix
  gslvector gslblock gslcomplex gslieeeutils gslsys gslerr utils gslrng
  gslsort)
