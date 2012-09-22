add_library(gsleigen
  jacobi.c symm.c symmv.c nonsymm.c nonsymmv.c herm.c hermv.c gensymm.c
  gensymmv.c genherm.c genhermv.c gen.c genv.c sort.c francis.c schur.c)

add_executable(gsleigen-test test.c)
target_link_libraries(gsleigen-test
  gsleigen gsltest gsllinalg gslpermutation gslblas gslcblas gslmatrix
  gslvector gslblock gslcomplex gslieeeutils gslsys gslerr utils gslrng
  gslsort)
if (NOT WIN32)
  target_link_libraries(gsleigen-test m)
endif ()
