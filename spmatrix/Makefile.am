noinst_LTLIBRARIES = libgslspmatrix.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_spmatrix.h gsl_spmatrix_char.h gsl_spmatrix_double.h gsl_spmatrix_float.h gsl_spmatrix_int.h gsl_spmatrix_long_double.h gsl_spmatrix_long.h gsl_spmatrix_short.h gsl_spmatrix_uchar.h gsl_spmatrix_uint.h gsl_spmatrix_ulong.h gsl_spmatrix_ushort.h gsl_spmatrix_complex_float.h gsl_spmatrix_complex_double.h gsl_spmatrix_complex_long_double.h

libgslspmatrix_la_SOURCES = compress.c copy.c file.c getset.c init.c minmax.c oper.c prop.c util.c swap.c

AM_CPPFLAGS = -I$(top_srcdir)

noinst_HEADERS = compress_source.c copy_source.c file_source.c getset_source.c getset_complex_source.c init_source.c minmax_source.c oper_source.c oper_complex_source.c prop_source.c swap_source.c test_source.c test_complex_source.c

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslspmatrix.la ../bst/libgslbst.la ../test/libgsltest.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../block/libgslblock.la  ../sys/libgslsys.la ../err/libgslerr.la ../utils/libutils.la ../rng/libgslrng.la
