add_library(gslieeeutils print.c make_rep.c gsl_ieee_utils.h env.c fp.c read.c)

add_executable(gslieeeutils-test test.c)
target_link_libraries(gslieeeutils-test
  gslieeeutils gslerr gsltest gslsys utils ${CMAKE_REQUIRED_LIBRARIES})
