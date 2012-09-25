# CMake build script for the GNU Scientific Library.

cmake_minimum_required(VERSION 2.6)
project(GSL)

set(PACKAGE "gsl")
set(PACKAGE_NAME ${PACKAGE})
set(PACKAGE_TARNAME ${PACKAGE})
set(PACKAGE_BUGREPORT "")
set(PACKAGE_URL "")
set(LT_OBJDIR ".libs/")

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

if (NOT (VERSION MATCHES "\\+"))
  set(RELEASED /**/)
endif ()

set(PACKAGE_STRING "${PACKAGE} ${VERSION}")

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

check_include_files(stdio.h STDC_HEADERS)

# Check for IEEE arithmetic interface type.
if (CMAKE_SYSTEM_NAME MATCHES Linux)
  if (CMAKE_SYSTEM_PROCESSOR MATCHES sparc)
    set(HAVE_GNUSPARC_IEEE_INTERFACE 1)
  elseif (CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
    set(HAVE_GNUPPC_IEEE_INTERFACE 1)
  elseif (CMAKE_SYSTEM_PROCESSOR MATCHES 86)
    set(HAVE_GNUX86_IEEE_INTERFACE 1)
  endif ()
elseif (CMAKE_SYSTEM_NAME MATCHES SunOS)
  set(HAVE_SUNOS4_IEEE_INTERFACE 1)
elseif (CMAKE_SYSTEM_NAME MATCHES Solaris)
  set(HAVE_SOLARIS_IEEE_INTERFACE 1)
elseif (CMAKE_SYSTEM_NAME MATCHES hpux)
  set(HAVE_HPUX_IEEE_INTERFACE 1)
elseif (CMAKE_SYSTEM_NAME MATCHES Darwin)
  if (CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
    set(HAVE_DARWIN_IEEE_INTERFACE 1)
  elseif (CMAKE_SYSTEM_PROCESSOR MATCHES 86)
    set(HAVE_DARWIN86_IEEE_INTERFACE 1)
  endif ()
elseif (CMAKE_SYSTEM_NAME MATCHES NetBSD)
  set(HAVE_NETBSD_IEEE_INTERFACE 1)
elseif (CMAKE_SYSTEM_NAME MATCHES OpenBSD)
  set(HAVE_OPENBSD_IEEE_INTERFACE 1)
elseif (CMAKE_SYSTEM_NAME MATCHES FreeBSD)
  set(HAVE_FREEBSD_IEEE_INTERFACE 1)
endif ()

# Check for FPU_SETCW.
if (HAVE_GNUX86_IEEE_INTERFACE)
  check_c_source_compiles("
    #include <fpu_control.h>
    #ifndef _FPU_SETCW
    #include <i386/fpu_control.h>
    #define _FPU_SETCW(cw) __setfpucw(cw)
    #endif
    int main() { unsigned short mode = 0 ; _FPU_SETCW(mode); }"
    HAVE_FPU_SETCW)
  if (NOT HAVE_FPU_SETCW)
    set(HAVE_GNUX86_IEEE_INTERFACE FALSE)
  endif ()
endif ()

# Check for SSE extensions.
if (HAVE_GNUX86_IEEE_INTERFACE)
  check_c_source_compiles("
    #include <stdlib.h>
    #define _FPU_SETMXCSR(cw) asm volatile (\"ldmxcsr %0\" : : \"m\" (*&cw))
    int main() { unsigned int mode = 0x1f80 ; _FPU_SETMXCSR(mode); exit(0); }"
    HAVE_FPU_X86_SSE)
endif ()

include(CheckLibraryExists)
check_library_exists(m cos "" HAVE_LIBM)
if (HAVE_LIBM)
  set(CMAKE_REQUIRED_LIBRARIES m)
endif ()

# Check IEEE comparisons, whether "x != x" is true for NaNs.
set(TEST_SRC ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.c)
file(WRITE ${TEST_SRC} "
  #include <math.h>
  int main (void)
  {
     int status; double inf, nan;
     inf = exp(1.0e10);
     nan = inf / inf ;
     status = (nan == nan);
     exit (status);
  }")
try_run(RUN_RESULT COMPILE_RESULT ${CMAKE_BINARY_DIR} ${TEST_SRC}
      CMAKE_FLAGS -DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES})
if (RUN_RESULT EQUAL 0)
  set(HAVE_IEEE_COMPARISONS 1)
endif ()

# Check for IEEE denormalized arithmetic.
file(WRITE ${TEST_SRC} "
  #include <math.h>
  int main (void)
  {
     int i, status;
     volatile double z = 1e-308;
     for (i = 0; i < 5; i++) { z = z / 10.0 ; };
     for (i = 0; i < 5; i++) { z = z * 10.0 ; };
     status = (z == 0.0);
     exit (status);
  }")
try_run(RUN_RESULT COMPILE_RESULT ${CMAKE_BINARY_DIR} ${TEST_SRC}
      CMAKE_FLAGS -DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES})
if (RUN_RESULT EQUAL 0)
  set(HAVE_IEEE_DENORMALS 1)
endif ()

# Check for long double stdio.
file(WRITE ${TEST_SRC} "
  #include <stdlib.h>
  #include <stdio.h>
  int main (void)
  {
    const char * s = "5678.25"; long double x = 1.234 ;
    fprintf(stderr,"%Lg\n",x) ;
    sscanf(s, "%Lg", &x);
    if (x == 5678.25) {exit (0);} else {exit(1); }
  }")
try_run(RUN_RESULT COMPILE_RESULT ${CMAKE_BINARY_DIR} ${TEST_SRC})
if (RUN_RESULT EQUAL 0)
  set(HAVE_PRINTF_LONGDOUBLE 1)
endif ()

if (NOT CMAKE_COMPILER_IS_GNUCC)
  file(WRITE ${TEST_SRC} "
    #include <limits.h>
    int main (void) { return CHAR_MIN == 0; }")
  try_run(RUN_RESULT COMPILE_RESULT ${CMAKE_BINARY_DIR} ${TEST_SRC})
  if (RUN_RESULT EQUAL 1)
    set(__CHAR_UNSIGNED__ 1)
  endif ()
endif ()

file(READ ${GSL_SOURCE_DIR}/config.h.in CONFIG)

# HAVE_DOPRNT is not defined because all current systems have vprintf
# (http://bit.ly/ROIckb).
foreach (var inline HAVE_INLINE HAVE_C99_INLINE HIDE_INLINE_STATIC
             size_t volatile
             RELEASED
             ${HEADER_VARS} STDC_HEADERS
             HAVE_DOPRNT
             HAVE_GNUSPARC_IEEE_INTERFACE
             HAVE_GNUM68K_IEEE_INTERFACE
             HAVE_GNUPPC_IEEE_INTERFACE
             HAVE_GNUX86_IEEE_INTERFACE
             HAVE_SUNOS4_IEEE_INTERFACE
             HAVE_SOLARIS_IEEE_INTERFACE
             HAVE_HPUX11_IEEE_INTERFACE
             HAVE_HPUX_IEEE_INTERFACE
             HAVE_TRU64_IEEE_INTERFACE
             HAVE_IRIX_IEEE_INTERFACE
             HAVE_AIX_IEEE_INTERFACE
             HAVE_FREEBSD_IEEE_INTERFACE
             HAVE_OS2EMX_IEEE_INTERFACE
             HAVE_NETBSD_IEEE_INTERFACE
             HAVE_OPENBSD_IEEE_INTERFACE
             HAVE_DARWIN_IEEE_INTERFACE
             HAVE_DARWIN86_IEEE_INTERFACE
             HAVE_FPU_X86_SSE
             HAVE_IEEE_COMPARISONS
             HAVE_IEEE_DENORMALS
             HAVE_PRINTF_LONGDOUBLE
             __CHAR_UNSIGNED__)
  if (${var})
    string(REGEX REPLACE "#( *)undef ${var}"
      "#\\1define ${var} ${${var}}" CONFIG "${CONFIG}")
  else ()
    string(REGEX REPLACE "#( *)undef ${var}"
      "/* #\\1undef ${var} */" CONFIG "${CONFIG}")
  endif ()
endforeach ()

# Remember to put a definition in config.h.in for each of these.
include(CheckSymbolExists)
check_symbol_exists(EXIT_SUCCESS "stdlib.h" HAVE_EXIT_SUCCESS)
check_symbol_exists(EXIT_FAILURE "stdlib.h" HAVE_EXIT_FAILURE)
if (HAVE_EXIT_SUCCESS AND HAVE_EXIT_FAILURE)
  set(HAVE_EXIT_SUCCESS_AND_FAILURE 1)
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

check_symbol_exists(memcpy string.h HAVE_MEMCPY)
check_symbol_exists(memmove string.h HAVE_MEMMOVE)
check_symbol_exists(strdup string.h HAVE_STRDUP)
check_symbol_exists(strtol stdlib.h HAVE_STRTOL)
check_symbol_exists(strtoul stdlib.h HAVE_STRTOUL)
check_symbol_exists(vprintf stdio.h HAVE_VPRINTF)

foreach (var PACKAGE PACKAGE_NAME PACKAGE_TARNAME VERSION PACKAGE_VERSION
             PACKAGE_BUGREPORT PACKAGE_URL PACKAGE_STRING LT_OBJDIR
             HAVE_LIBM
             HAVE_EXIT_SUCCESS_AND_FAILURE
             HAVE_DECL_FEENABLEEXCEPT HAVE_DECL_FESETTRAPENABLE
             HAVE_DECL_HYPOT HAVE_DECL_EXPM1 HAVE_DECL_ACOSH
             HAVE_DECL_ASINH HAVE_DECL_ATANH HAVE_DECL_LDEXP
             HAVE_DECL_FREXP HAVE_DECL_FPRND_T HAVE_DECL_ISINF
             HAVE_DECL_ISFINITE HAVE_DECL_FINITE HAVE_DECL_ISNAN
             HAVE_DECL_LOG1P
             HAVE_EXTENDED_PRECISION_REGISTERS
             HAVE_MEMCPY HAVE_MEMMOVE HAVE_STRDUP HAVE_STRTOL HAVE_STRTOUL
             HAVE_VPRINTF)
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
