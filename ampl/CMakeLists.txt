# CMake build script for the GSL wrapper for AMPL.
if(MSVC_RUNTIME_DYNAMIC)
  set(BUILD_DYNRT_LIBS)
  set(asllib ${asl-dynrt})
else()
  set(asllib ${asl})
endif()
add_subdirectory(thirdparty/asl)

set(AMPLGSL_VERSION 20211014)
# Install location
set(AMPL_LIBRARY_DIR gsl)

if (MSVC)
  # Disable useless MSVC warnings suggesting nonportable "secure" alternatives.
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
endif ()

# Package name definition
if(LINUX)
  set(SYSNAME "linux-intel")
elseif(APPLE)
  set(SYSNAME "macos")
elseif(WIN32)
  set(SYSNAME "mswin")
endif()
set(PACKAGENAME "amplgsl.${SYSNAME}${ARCH}.${AMPLGSL_VERSION}")

# Macro to build an AMPL library (links with ASL and has dll extension)
macro(add_ampl_library name)
  cmake_parse_arguments(add_ampl_library PRIVATE "" "" ${ARGN})
  add_library(${name} SHARED ${add_ampl_library_UNPARSED_ARGUMENTS})
  set_target_properties(${name} PROPERTIES PREFIX "")
  set_target_properties(${name} PROPERTIES SUFFIX ".dll")
  
  target_link_libraries(${name} ${asllib})
  if(WIN32)
    set(toInstall "RUNTIME")
  else()
    set(toInstall "LIBRARY")
  endif()
  # We don't want to import libraries to be installed.
   install(TARGETS ${name} ${toInstall}
   DESTINATION "."
    COMPONENT ${PACKAGENAME})
endmacro()

# Output GSL version number
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gsl-version.txt "AMPLGSL version ${AMPLGSL_VERSION}\nGSL version ${VERSION}")

# Definition of AMPL-gsl interface
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_ampl_library(amplgsl ${SRCDIR}/amplgsl.cc)
target_link_libraries(amplgsl gsl gslcblas asl)
target_include_directories(amplgsl PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/asl/src/solvers)

# Used to generate the ampl function definitions
add_executable(gsl-info ${SRCDIR}/gsl-info.cc)
target_compile_definitions(gsl-info PRIVATE LIBDATE=${AMPLGSL_VERSION})
target_link_libraries(gsl-info amplgsl)
add_custom_command(OUTPUT gsl.ampl COMMAND gsl-info DEPENDS amplgsl gsl-info)
add_custom_target(gsl-ampl ALL DEPENDS gsl.ampl)
# Install AMPL script
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gsl.ampl ${CMAKE_CURRENT_BINARY_DIR}/gsl-version.txt
    DESTINATION "."
    COMPONENT ${PACKAGENAME})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/CHANGES.amplgsl.md
    DESTINATION "docs"
    COMPONENT ${PACKAGENAME})

 # Helper library from mp
add_subdirectory(thirdparty/test-support)

# Test target
add_executable(amplgsl-test test/gsl-test.cc)
if (MINGW)
   set_target_properties(amplgsl-test PROPERTIES
     LINK_FLAGS "-static-libgcc -static-libstdc++")
endif()
target_link_libraries(amplgsl-test amplgsl test-support)
if(MSVC)
    target_compile_options(amplgsl-test PRIVATE /bigobj)
endif()
target_compile_definitions(amplgsl-test PRIVATE
    AMPLGSL_DLL_NAME="$<TARGET_FILE:amplgsl>")
add_test(NAME amplgsl-test COMMAND $<TARGET_FILE:amplgsl-test>)


# Documentation
add_subdirectory(doc)

######  Package creation settings ######
# Packages definitions happens in the setSolverProperties function
# in file SolversListFunctions.cmake
if(MSVC)
    set(CPACK_GENERATOR ZIP)
else()
    set(CPACK_GENERATOR TGZ)
endif()
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY OFF)
set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)
set(CPACK_PACKAGE_VERSION 1.0})
set(CPACK_PACKAGE_FILE_NAME amplgsl)
include(CPack)
