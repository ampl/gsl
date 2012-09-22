add_library(gslcblas
  sasum.c saxpy.c scasum.c scnrm2.c scopy.c sdot.c sdsdot.c sgbmv.c sgemm.c
  sgemv.c sger.c snrm2.c srot.c srotg.c srotm.c srotmg.c ssbmv.c sscal.c
  sspmv.c sspr.c sspr2.c sswap.c ssymm.c ssymv.c ssyr.c ssyr2.c ssyr2k.c
  ssyrk.c stbmv.c stbsv.c stpmv.c stpsv.c strmm.c strmv.c strsm.c strsv.c
  dasum.c daxpy.c dcopy.c ddot.c dgbmv.c dgemm.c dgemv.c dger.c dnrm2.c
  drot.c drotg.c drotm.c drotmg.c dsbmv.c dscal.c dsdot.c dspmv.c dspr.c
  dspr2.c dswap.c dsymm.c dsymv.c dsyr.c dsyr2.c dsyr2k.c dsyrk.c dtbmv.c
  dtbsv.c dtpmv.c dtpsv.c dtrmm.c dtrmv.c dtrsm.c dtrsv.c dzasum.c dznrm2.c
  caxpy.c ccopy.c cdotc_sub.c cdotu_sub.c cgbmv.c cgemm.c cgemv.c cgerc.c
  cgeru.c chbmv.c chemm.c chemv.c cher.c cher2.c cher2k.c cherk.c chpmv.c
  chpr.c chpr2.c cscal.c csscal.c cswap.c csymm.c csyr2k.c csyrk.c ctbmv.c
  ctbsv.c ctpmv.c ctpsv.c ctrmm.c ctrmv.c ctrsm.c ctrsv.c zaxpy.c zcopy.c
  zdotc_sub.c zdotu_sub.c zdscal.c zgbmv.c zgemm.c zgemv.c zgerc.c zgeru.c
  zhbmv.c zhemm.c zhemv.c zher.c zher2.c zher2k.c zherk.c zhpmv.c zhpr.c
  zhpr2.c zscal.c zswap.c zsymm.c zsyr2k.c zsyrk.c ztbmv.c ztbsv.c ztpmv.c
  ztpsv.c ztrmm.c ztrmv.c ztrsm.c ztrsv.c icamax.c idamax.c isamax.c izamax.c
  xerbla.c)

add_executable(gslcblas-test
  test.c test_amax.c test_asum.c test_axpy.c test_copy.c test_dot.c
  test_gbmv.c test_gemm.c test_gemv.c test_ger.c test_hbmv.c test_hemm.c
  test_hemv.c test_her.c test_her2.c test_her2k.c test_herk.c test_hpmv.c
  test_hpr.c test_hpr2.c test_nrm2.c test_rot.c test_rotg.c test_rotm.c
  test_rotmg.c test_sbmv.c test_scal.c test_spmv.c test_spr.c test_spr2.c
  test_swap.c test_symm.c test_symv.c test_syr.c test_syr2.c test_syr2k.c
  test_syrk.c test_tbmv.c test_tbsv.c test_tpmv.c test_tpsv.c test_trmm.c
  test_trmv.c test_trsm.c test_trsv.c)
target_link_libraries(gslcblas-test
  gslcblas gslieeeutils gslerr gsltest gslsys)
if (NOT WIN32)
  target_link_libraries(gslcblas-test m)
endif ()
