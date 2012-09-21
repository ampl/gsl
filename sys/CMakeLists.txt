add_library(gslsys
  minmax.c prec.c hypot.c log1p.c expm1.c coerce.c
  invhyp.c pow_int.c infnan.c fdiv.c fcmp.c ldfrexp.c)

add_executable(gslsys-test test.c)
target_link_libraries(gslsys-test
  gslsys gslieeeutils gslerr gsltest gslsys utils)
if (NOT WIN32)
  target_link_libraries(gslsys-test m)
endif ()
