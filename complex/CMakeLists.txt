add_library(gslcomplex math.c inline.c)

add_executable(gslcomplex-test
  test.c results.h results1.h results2.h results_real.h results_zreal.h)
target_link_libraries(gslcomplex-test
  gslcomplex gslerr gslieeeutils gslerr gsltest gslsys utils
  ${CMAKE_REQUIRED_LIBRARIES})
