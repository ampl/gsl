add_library(gslerr error.c stream.c message.c strerror.c)

add_executable(gslerr-test test.c)
target_link_libraries(gslerr-test gslerr gsltest gslsys utils)
