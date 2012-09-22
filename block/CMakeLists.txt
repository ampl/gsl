add_library(gslblock init.c file.c block.c)

add_executable(gslblock-test test.c)
target_link_libraries(gslblock-test
  gslblock gslieeeutils gslerr gsltest gslsys utils)
if (NOT WIN32)
  target_link_libraries(gslblock-test m)
endif ()
