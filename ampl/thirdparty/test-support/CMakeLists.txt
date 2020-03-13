  # Helper library from mp
  find_package(Threads)
  include(CheckCXXCompilerFlag)
  check_cxx_source_compiles("
    #include <initializer_list>
    int main() {}" MP_INITIALIZER_LISTS)
  if (NOT MP_VARIADIC_TEMPLATES OR NOT MP_INITIALIZER_LISTS)
    add_definitions(-DGTEST_LANG_CXX11=0)
  endif ()
  set(MP_TEST_SOURCES 
      ${CMAKE_CURRENT_SOURCE_DIR}/function.h
      ${CMAKE_CURRENT_SOURCE_DIR}/util.h
      ${CMAKE_CURRENT_SOURCE_DIR}/function.cc
      ${CMAKE_CURRENT_SOURCE_DIR}/format.cc
      ${CMAKE_CURRENT_SOURCE_DIR}/util.cc
      ${CMAKE_CURRENT_SOURCE_DIR}/test-main.cc
      ${CMAKE_CURRENT_SOURCE_DIR}/../gmock/gmock-gtest-all.cc)
  add_library(test-support STATIC ${MP_TEST_SOURCES})
  
  target_link_libraries(test-support PRIVATE ${CMAKE_THREAD_LIBS_INIT} asl2)
  target_include_directories(test-support PUBLIC 
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/..)
  # Workaround GTest bug https://github.com/google/googletest/issues/705.
  check_cxx_compiler_flag(
    -fno-delete-null-pointer-checks HAVE_FNO_DELETE_NULL_POINTER_CHECKS)
  if (HAVE_FNO_DELETE_NULL_POINTER_CHECKS)
    target_compile_options(test-support PUBLIC -fno-delete-null-pointer-checks)
  endif ()
  # Workaround a bug in implementation of variadic templates in MSVC11.
  if (MSVC)
    target_compile_definitions(test-support PUBLIC _VARIADIC_MAX=10)
  endif ()
