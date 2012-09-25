add_library(gslsort
  sort.c sortind.c sortvec.c sortvecind.c subset.c subsetind.c)

add_executable(gslsort-test test.c)
target_link_libraries(gslsort-test
  gslsort gslpermutation gslvector gslblock gslieeeutils gslerr gsltest gslsys
  ${CMAKE_REQUIRED_LIBRARIES})
