add_library(gslpermutation
  init.c file.c permutation.c permute.c canonical.c inline.c)

add_executable(gslpermutation-test test.c)
target_link_libraries(gslpermutation-test
  gslpermutation gslvector gslblock gslieeeutils gslerr gsltest gslsys utils
  ${CMAKE_REQUIRED_LIBRARIES})
