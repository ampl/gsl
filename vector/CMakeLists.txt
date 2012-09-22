add_library(gslvector
  init.c file.c vector.c copy.c swap.c prop.c minmax.c oper.c
  reim.c subvector.c view.c view.h)

add_executable(gslvector-test test.c)
target_link_libraries(gslvector-test
  gslvector gslblock gslieeeutils gslerr gsltest gslsys utils)
if (NOT WIN32)
  target_link_libraries(gslvector-test m)
endif ()

add_executable(gslvector-test-static test_static.c)
target_link_libraries(gslvector-test-static
  gslvector gslblock gslieeeutils gslerr gsltest gslsys utils)
if (NOT WIN32)
  target_link_libraries(gslvector-test-static m)
endif ()
