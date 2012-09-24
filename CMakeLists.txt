# CMake build script for the GNU Scientific Library.

cmake_minimum_required(VERSION 2.6)
project(GSL)

set(PACKAGE "gsl")
set(PACKAGE_NAME ${PACKAGE})
set(PACKAGE_TARNAME ${PACKAGE})
set(PACKAGE_BUGREPORT "")
set(PACKAGE_URL "")

# Get version numbers from configure.ac.
file(STRINGS configure.ac CONFIG REGEX
  "^(AC_INIT|(GSL|CBLAS)_(CURRENT|REVISION|AGE)=)")
foreach (line ${CONFIG})
  if (line MATCHES AC_INIT)
    string(REGEX REPLACE "AC_INIT\\(\\[gsl\\],\\[(.*)\\]\\)" "\\1"
      VERSION ${line})
  else ()
    string(REGEX REPLACE "(.*)=.*" "\\1" name ${line})
    string(REGEX REPLACE ".*=(.*)" "\\1" value ${line})
    set(${name} ${value})
  endif ()
endforeach ()
set(PACKAGE_VERSION ${VERSION})

include(CheckCSourceCompiles)

foreach (keyword inline __inline__ __inline)
  check_c_source_compiles("
    static ${keyword} void foo() { return 0; }
    int main() {}" C_HAS_${keyword})
  if (C_HAS_${keyword})
    set(C_INLINE ${keyword})
    break ()
  endif ()
endforeach ()
add_definitions("-Dinline=${C_INLINE}")
if (C_INLINE)
  # Check for GNU-style extern inline.
  check_c_source_compiles("
    extern ${C_INLINE} double foo(double x);
    extern ${C_INLINE} double foo(double x) { return x + 1.0; }
    double foo(double x) { return x + 1.0; }
    int main() { foo(1.0); }" C_EXTERN_INLINE)
  if (C_EXTERN_INLINE)
    set(HAVE_INLINE 1)
  else ()
    # Check for C99-style inline.
    check_c_source_compiles("
      extern inline void* foo() { foo(); return &foo; }
      int main() { return foo() != 0; }" C_C99INLINE)
    if (C_C99INLINE)
      set(HAVE_INLINE 1)
      set(HAVE_C99_INLINE 1)
    endif ()
  endif ()
endif ()
if (C_INLINE AND NOT C_HAS_inline)
  set(inline ${C_INLINE})
endif ()

include(CheckIncludeFiles)

# Checks for header files.
macro(check_headers headers)
  foreach (header ${headers})
    string(TOUPPER HAVE_${header} var)
    string(REGEX REPLACE "\\.|/" "_" var ${var})
    check_include_files(${header} ${var})
    set(HEADER_VARS ${HEADER_VARS} ${var})
  endforeach ()
endmacro()

set(HEADERS
  ieeefp.h dlfcn.h inttypes.h memory.h stdint.h stdlib.h
  strings.h string.h sys/stat.h sys/types.h unistd.h)
check_headers("${HEADERS}")

file(READ ${GSL_SOURCE_DIR}/config.h.in CONFIG)

# HAVE_DOPRNT is not defined because all current systems have vprintf
# (http://bit.ly/ROIckb).
foreach (var inline HAVE_INLINE HAVE_C99_INLINE HAVE_DOPRNT ${HEADER_VARS})
  if (${var})
    string(REPLACE "#undef ${var}"
      "#define ${var} ${${var}}" CONFIG "${CONFIG}")
  else ()
    string(REPLACE "#undef ${var}"
      "/* #undef ${var} */" CONFIG "${CONFIG}")
  endif ()
endforeach ()

# Remember to put a definition in config.h.in for each of these.
include(CheckSymbolExists)
check_symbol_exists(EXIT_SUCCESS "stdlib.h" HAVE_EXIT_SUCCESS)
check_symbol_exists(EXIT_FAILURE "stdlib.h" HAVE_EXIT_FAILURE)
if (HAVE_EXIT_SUCCESS AND HAVE_EXIT_FAILURE)
  set(HAVE_EXIT_SUCCESS_AND_FAILURE 1)
endif ()
if (NOT WIN32)
  set(CMAKE_REQUIRED_LIBRARIES m)
endif ()
set(CMAKE_REQUIRED_DEFINITIONS "-D_GNU_SOURCE=1")
check_symbol_exists(feenableexcept fenv.h HAVE_DECL_FEENABLEEXCEPT)
check_symbol_exists(fesettrapenable fenv.h HAVE_DECL_FESETTRAPENABLE)
set(CMAKE_REQUIRED_DEFINITIONS "")
check_symbol_exists(hypot math.h HAVE_DECL_HYPOT)
check_symbol_exists(expm1 math.h HAVE_DECL_EXPM1)
check_symbol_exists(acosh math.h HAVE_DECL_ACOSH)
check_symbol_exists(asinh math.h HAVE_DECL_ASINH)
check_symbol_exists(atanh math.h HAVE_DECL_ATANH)
check_symbol_exists(ldexp math.h HAVE_DECL_LDEXP)
check_symbol_exists(frexp math.h HAVE_DECL_FREXP)
check_symbol_exists(fprnd_t float.h HAVE_DECL_FPRND_T)
check_symbol_exists(isinf math.h HAVE_DECL_ISINF)
check_symbol_exists(isfinite math.h HAVE_DECL_ISFINITE)
if (HAVE_IEEEFP_H)
  set(IEEEFP_H ieeefp.h)
endif ()
check_symbol_exists(finite math.h;${IEEEFP_H} HAVE_DECL_FINITE)
check_symbol_exists(isnan math.h HAVE_DECL_ISNAN)

# OpenBSD has a broken implementation of log1p.
if (${CMAKE_SYSTEM_NAME} MATCHES "OpenBSD")
  message("avoiding OpenBSD system log1p - using gsl version")
else ()
  check_symbol_exists(log1p math.h HAVE_DECL_LOG1P)
endif ()

# Check for extended floating point registers.
if (NOT (CMAKE_SYSTEM_PROCESSOR MATCHES "^(sparc|powerpc|hppa|alpha)"))
  set(HAVE_EXTENDED_PRECISION_REGISTERS 1)
endif ()

foreach (var PACKAGE PACKAGE_NAME PACKAGE_TARNAME VERSION PACKAGE_VERSION
             PACKAGE_BUGREPORT PACKAGE_URL
             HAVE_EXIT_SUCCESS_AND_FAILURE
             HAVE_DECL_FEENABLEEXCEPT HAVE_DECL_FESETTRAPENABLE
             HAVE_DECL_HYPOT HAVE_DECL_EXPM1 HAVE_DECL_ACOSH
             HAVE_DECL_ASINH HAVE_DECL_ATANH HAVE_DECL_LDEXP
             HAVE_DECL_FREXP HAVE_DECL_FPRND_T HAVE_DECL_ISINF
             HAVE_DECL_ISFINITE HAVE_DECL_FINITE HAVE_DECL_ISNAN
             HAVE_DECL_LOG1P
             HAVE_EXTENDED_PRECISION_REGISTERS)
  if (NOT (var MATCHES ^HAVE))
    set(value \"${${var}}\")
  elseif (${var})
    set(value 1)
  else ()
    set(value 0)
  endif ()
  string(REGEX REPLACE "#undef ${var}([\r\n])"
    "#define ${var} ${value}\\1" CONFIG "${CONFIG}")
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
