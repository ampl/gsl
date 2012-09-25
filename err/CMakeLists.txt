add_library(gslerr error.c stream.c message.c strerror.c)

add_gsl_test(gslerr-test test.c LIBS gslerr gsltest gslsys utils)
