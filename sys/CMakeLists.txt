add_library(gslsys
  minmax.c prec.c hypot.c log1p.c expm1.c coerce.c
  invhyp.c pow_int.c infnan.c fdiv.c fcmp.c ldfrexp.c)

add_gsl_test(gslsys-test test.c
  LIBS gslsys gslieeeutils gslerr gsltest gslsys utils)
