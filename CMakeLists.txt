# CMake build script for the GNU Scientific Library.

cmake_minimum_required(VERSION 2.6)
project(GSL)

# Get version numbers from configure.ac.
file(STRINGS configure.ac CONFIG REGEX
  "^(GSL|CBLAS)_(CURRENT|REVISION|AGE)=.*")
foreach (line ${CONFIG})
  string(REGEX REPLACE "(.*)=.*" "\\1" name ${line})
  string(REGEX REPLACE ".*=(.*)" "\\1" value ${line})
  set(${name} ${value})
endforeach ()

include(CheckSymbolExists)
include(CheckCXXSourceCompiles)

file(READ ${GSL_SOURCE_DIR}/config.h.in CONFIG)

# Checks for header files.
include(CheckIncludeFiles)
check_include_files(ieeefp.h HAVE_IEEEFP_H)
foreach (var HAVE_IEEEFP_H)
  if (${var})
    string(REPLACE "#undef ${var}"
      "#define ${var} ${value}" CONFIG "${CONFIG}")
  else ()
    string(REPLACE "#undef ${var}"
      "/* #undef ${var} */" CONFIG "${CONFIG}")
  endif ()
endforeach ()

# Remember to put a definition in config.h.in for each of these.
if (NOT WIN32)
  set(CMAKE_REQUIRED_LIBRARIES m)
endif ()
set(CMAKE_REQUIRED_DEFINITIONS "-D_GNU_SOURCE=1")
check_symbol_exists(feenableexcept "fenv.h" HAVE_DECL_FEENABLEEXCEPT)
check_symbol_exists(fesettrapenable "fenv.h" HAVE_DECL_FESETTRAPENABLE)
set(CMAKE_REQUIRED_DEFINITIONS "")
check_symbol_exists(hypot "math.h" HAVE_DECL_HYPOT)
check_symbol_exists(expm1 "math.h" HAVE_DECL_EXPM1)
check_symbol_exists(acosh "math.h" HAVE_DECL_ACOSH)
check_symbol_exists(asinh "math.h" HAVE_DECL_ASINH)
check_symbol_exists(atanh "math.h" HAVE_DECL_ATANH)
check_symbol_exists(ldexp "math.h" HAVE_DECL_LDEXP)
check_symbol_exists(frexp "math.h" HAVE_DECL_FREXP)
check_symbol_exists(fprnd_t "float.h" HAVE_DECL_FPRND_T)
check_symbol_exists(isinf "math.h" HAVE_DECL_ISINF)
check_symbol_exists(isfinite "math.h" HAVE_DECL_ISFINITE)
if (HAVE_IEEEFP_H)
  set(IEEEFP_H ieeefp.h)
endif ()
check_symbol_exists(finite "math.h;${IEEEFP_H}" HAVE_DECL_FINITE)
check_symbol_exists(isnan "math.h" HAVE_DECL_ISNAN)

# OpenBSD has a broken implementation of log1p.
if (${CMAKE_SYSTEM_NAME} MATCHES "OpenBSD")
  message("avoiding OpenBSD system log1p - using gsl version")
else ()
  check_symbol_exists(log1p "math.h" HAVE_DECL_LOG1P)
endif ()

foreach (var HAVE_DECL_FEENABLEEXCEPT HAVE_DECL_FESETTRAPENABLE
             HAVE_DECL_HYPOT HAVE_DECL_EXPM1 HAVE_DECL_ACOSH
             HAVE_DECL_ASINH HAVE_DECL_ATANH HAVE_DECL_LDEXP
             HAVE_DECL_FREXP HAVE_DECL_FPRND_T HAVE_DECL_ISINF
             HAVE_DECL_ISFINITE HAVE_DECL_FINITE HAVE_DECL_ISNAN
             HAVE_DECL_LOG1P)
  set(value 0)
  if (${var})
    set(value 1)
  endif ()
  string(REPLACE "#undef ${var}"
    "#define ${var} ${value}" CONFIG "${CONFIG}")
endforeach ()

file(WRITE ${GSL_BINARY_DIR}/config.h
"/* config.h.  Generated from config.h.in by configure.  */
${CONFIG}")

include_directories(${GSL_BINARY_DIR} ${GSL_SOURCE_DIR})

add_subdirectory(gsl)
add_subdirectory(utils)
add_subdirectory(sys)
add_subdirectory(test)
add_subdirectory(err)
add_subdirectory(complex)
add_subdirectory(block)
add_subdirectory(vector)
add_subdirectory(matrix)
add_subdirectory(permutation)
add_subdirectory(sort)
add_subdirectory(ieee-utils)
add_subdirectory(cblas)
add_subdirectory(blas)
add_subdirectory(linalg)
add_subdirectory(eigen)
add_subdirectory(specfunc)
add_subdirectory(rng)
