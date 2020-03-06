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
