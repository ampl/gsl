add_library(gslmatrix
  init.c matrix.c file.c rowcol.c swap.c copy.c minmax.c prop.c
  oper.c getset.c view.c submatrix.c view.h)

add_executable(gslmatrix-test test.c)
target_link_libraries(gslmatrix-test
  gslmatrix gslvector gslblock gslieeeutils gslerr gsltest gslsys utils)
if (NOT WIN32)
  target_link_libraries(gslmatrix-test m)
endif ()

add_executable(gslmatrix-test-static test_static.c)
target_link_libraries(gslmatrix-test-static
  gslmatrix gslvector gslblock gslieeeutils gslerr gsltest gslsys utils)
if (NOT WIN32)
  target_link_libraries(gslmatrix-test-static m)
endif ()
