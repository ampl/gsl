# CMake build script for the GSL wrapper for AMPL.
add_subdirectory(thirdparty/asl)
# Install location
set(AMPL_LIBRARY_DIR gsl)

# Macro to build an AMPL library (links with ASL and has dll extension)
macro(add_ampl_library name)
  cmake_parse_arguments(add_ampl_library PRIVATE "" "" ${ARGN})
  add_library(${name} SHARED ${add_ampl_library_UNPARSED_ARGUMENTS})
  set_target_properties(${name} PROPERTIES PREFIX "")
  set_target_properties(${name} PROPERTIES SUFFIX ".dll")
  target_link_libraries(${name} asl)
  # Specify RUNTIME DESTINATION and LIBRARY DESTINATION, but not
  # DESTINATION or ARCHIVE_DESTINATION because we don't want to import
  # libraries installed.
   install(TARGETS ${name} 
    RUNTIME DESTINATION ${AMPL_LIBRARY_DIR} 
    LIBRARY DESTINATION ${AMPL_LIBRARY_DIR})
endmacro()

# Definition of AMPL-gsl interface
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_ampl_library(amplgsl ${SRCDIR}/amplgsl.cc)
target_link_libraries(amplgsl gsl gslcblas asl)
target_include_directories(amplgsl PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/asl/src/solvers)

# Used to generate the ampl function definitions
add_executable(gsl-info ${SRCDIR}/gsl-info.cc)
target_link_libraries(gsl-info amplgsl)
add_custom_command(OUTPUT gsl.ampl COMMAND gsl-info DEPENDS amplgsl gsl-info)
add_custom_target(gsl-ampl ALL DEPENDS gsl.ampl)
# Install AMPL script
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gsl.ampl 
    DESTINATION ${AMPL_LIBRARY_DIR})

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
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/mptest/function.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/mptest/format.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/mptest/test-main.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/gmock/gmock-gtest-all.cc)
add_library(amplgsl-test-support STATIC ${MP_TEST_SOURCES})
target_include_directories(amplgsl-test-support PRIVATE 
   ${ASL1_INCLUDE_DIRS})
target_link_libraries(amplgsl-test-support PRIVATE ${CMAKE_THREAD_LIBS_INIT})
target_include_directories(amplgsl-test-support PUBLIC 
   ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/mptest
   ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
# Workaround GTest bug https://github.com/google/googletest/issues/705.
check_cxx_compiler_flag(
  -fno-delete-null-pointer-checks HAVE_FNO_DELETE_NULL_POINTER_CHECKS)
if (HAVE_FNO_DELETE_NULL_POINTER_CHECKS)
  target_compile_options(amplgsl-test-support PUBLIC -fno-delete-null-pointer-checks)
endif ()
# Workaround a bug in implementation of variadic templates in MSVC11.
if (MSVC)
  target_compile_definitions(amplgsl-test-support PUBLIC _VARIADIC_MAX=10)
endif ()

# Test target
add_executable(amplgsl-test test/gsl-test.cc)
if (MINGW)
   set_target_properties(amplgsl-test PROPERTIES
     LINK_FLAGS "-static-libgcc -static-libstdc++")
endif()
target_link_libraries(amplgsl-test amplgsl amplgsl-test-support)
if(MSVC)
    target_compile_options(amplgsl-test PRIVATE /bigobj)
endif()
target_compile_definitions(amplgsl-test PRIVATE
    AMPLGSL_DLL_NAME="$<TARGET_FILE:amplgsl>")
add_test(NAME amplgsl-test COMMAND $<TARGET_FILE:amplgsl-test>)
