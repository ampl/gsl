add_library(gsllinalg
  multiply.c exponential.c tridiag.c tridiag.h lu.c luc.c hh.c qr.c qrpt.c
  lq.c ptlq.c svd.c householder.c householdercomplex.c hessenberg.c hesstri.c
  cholesky.c choleskyc.c symmtd.c hermtd.c bidiag.c balance.c balancemat.c)

add_executable(gsllinalg-test test.c)
target_link_libraries(gsllinalg-test
  gsllinalg gslblas gslcblas gslpermutation gslmatrix gslvector
  gslblock gslcomplex gslieeeutils gslerr gsltest gslsys utils
  ${CMAKE_REQUIRED_LIBRARIES})
