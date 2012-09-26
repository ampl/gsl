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
file(STRINGS configure.ac LINES)
foreach (line ${LINES})
  if (line MATCHES "AC_INIT")
    string(REGEX REPLACE "AC_INIT[^,]*,\\[(.*)\\].*" "\\1" VERSION ${line})
    message("Got VERSION=${VERSION} from configure.ac")
  elseif (line MATCHES "(GSL|CBLAS)_(CURRENT|REVISION|AGE)=")
    string(REGEX REPLACE "(.*)=.*" "\\1" name ${line})
    string(REGEX REPLACE ".*=(.*)" "\\1" value ${line})
    message("Got ${name}=${value} from configure.ac")
    set(${name} ${value})
  endif ()
endforeach ()
set(PACKAGE_VERSION ${VERSION})
set(PACKAGE_STRING "${PACKAGE} ${VERSION}")

if (NOT (VERSION MATCHES "\\+"))
  # Defined if this is an official release.
  set(RELEASED /**/)
endif ()

include(CheckLibraryExists)
check_library_exists(m cos "" HAVE_LIBM)
if (HAVE_LIBM)
  set(CMAKE_REQUIRED_LIBRARIES m)
endif ()

include(CheckCSourceCompiles)

# Check for inline.
foreach (keyword inline __inline__ __inline)
  check_c_source_compiles("
    static ${keyword} void foo() { return 0; }
    int main() {}" C_HAS_${keyword})
  if (C_HAS_${keyword})
    set(C_INLINE ${keyword})
    break ()
  endif ()
endforeach ()
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

# Checks for header files.
include(CheckIncludeFiles)
foreach (header ieeefp.h dlfcn.h inttypes.h memory.h stdint.h stdlib.h
                strings.h string.h sys/stat.h sys/types.h unistd.h)
  string(TOUPPER HAVE_${header} var)
  string(REGEX REPLACE "\\.|/" "_" var ${var})
  check_include_files(${header} ${var})
endforeach ()
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
    set(HAVE_GNUX86_IEEE_INTERFACE 0)
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

# Compiles the source code, runs the program and sets ${VAR} to 1 if the
# return value is equal to ${RESULT}.
macro(check_run_result SRC RESULT VAR)
  set(SRC_FILE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.c)
  file(WRITE ${SRC_FILE} "${SRC}")
  try_run(RUN_RESULT COMPILE_RESULT ${CMAKE_BINARY_DIR} ${SRC_FILE}
          CMAKE_FLAGS -DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES})
  if (RUN_RESULT EQUAL ${RESULT})
    set(${VAR} 1)
  endif ()
endmacro()

# Check IEEE comparisons, whether "x != x" is true for NaNs.
check_run_result("
  #include <math.h>
  int main (void)
  {
     int status; double inf, nan;
     inf = exp(1.0e10);
     nan = inf / inf ;
     status = (nan == nan);
     exit (status);
  }" 0 HAVE_IEEE_COMPARISONS)

# Check for IEEE denormalized arithmetic.
check_run_result("
  #include <math.h>
  int main (void)
  {
     int i, status;
     volatile double z = 1e-308;
     for (i = 0; i < 5; i++) { z = z / 10.0 ; };
     for (i = 0; i < 5; i++) { z = z * 10.0 ; };
     status = (z == 0.0);
     exit (status);
  }" 0 HAVE_IEEE_DENORMALS)

# Check for long double stdio.
check_run_result("
  #include <stdlib.h>
  #include <stdio.h>
  int main (void)
  {
    const char * s = \"5678.25\"; long double x = 1.234 ;
    fprintf(stderr,\"%Lg\n\",x) ;
    sscanf(s, \"%Lg\", &x);
    if (x == 5678.25) {exit (0);} else {exit(1); }
  }" 0 HAVE_PRINTF_LONGDOUBLE)

if (NOT CMAKE_COMPILER_IS_GNUCC)
  check_run_result("
    #include <limits.h>
    int main (void) { return CHAR_MIN == 0; }" 1 __CHAR_UNSIGNED__)
endif ()

# Remember to put a definition in config.h.in for each of these.
include(CheckSymbolExists)
check_symbol_exists(EXIT_SUCCESS stdlib.h HAVE_EXIT_SUCCESS)
check_symbol_exists(EXIT_FAILURE stdlib.h HAVE_EXIT_FAILURE)
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
if (CMAKE_SYSTEM_NAME MATCHES OpenBSD)
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

# Process config.h using autoconf rules.
file(STRINGS ${GSL_SOURCE_DIR}/config.h.in CONFIG)
list(LENGTH CONFIG length)
math(EXPR length "${length} - 1")
foreach (i RANGE ${length})
  list(GET CONFIG ${i} line)
  if (line MATCHES "^# *undef ")
    string(REGEX REPLACE "# *undef (.*)" "\\1" var ${line})
    if (NOT DEFINED ${var} OR (var MATCHES "HAVE_.*_H" AND NOT ${var}))
      string(REGEX REPLACE "#( *)undef ${var}"
             "/* #\\1undef ${var} */" line "${line}")
    else ()
      if ("${${var}}" STREQUAL "/**/" OR "${var}" STREQUAL "inline")
        set(value ${${var}})
      elseif (NOT (var MATCHES ^HAVE OR ${var} EQUAL 0 OR ${var} EQUAL 1))
        set(value \"${${var}}\")
      elseif (${var})
        set(value 1)
      else ()
        set(value 0)
      endif ()
      string(REGEX REPLACE "#( *)undef ${var}"
             "#\\1define ${var} ${value}" line "${line}")
    endif ()
  endif ()
  set(CONFIG_OUT "${CONFIG_OUT}${line}\n")
endforeach ()
file(WRITE ${GSL_BINARY_DIR}/config.h
"/* config.h.  Generated from config.h.in by configure.  */
${CONFIG_OUT}")

include_directories(${GSL_BINARY_DIR} ${GSL_SOURCE_DIR})
add_definitions(-DHAVE_CONFIG_H)

enable_testing()

# Adds a GSL test. Usage:
#   add_gsl_test(<exename> <source> ...)
function(add_gsl_test exename)
  add_executable(${exename} ${ARGN})
  target_link_libraries(${exename} gsl gslcblas ${CMAKE_REQUIRED_LIBRARIES})
  add_test(${exename} ${exename})
endfunction()

file(GLOB GSL_HEADER_PATHS "${CMAKE_SOURCE_DIR}/gsl*.h"
  "${CMAKE_SOURCE_DIR}/*/gsl*.h")

foreach (path ${GSL_HEADER_PATHS})
  get_filename_component(name ${path} NAME)
  set(GSL_HEADERS ${GSL_HEADERS} ${name})
endforeach ()

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/copy-headers.cmake "
  file(MAKE_DIRECTORY gsl)
  foreach (path ${GSL_HEADER_PATHS})
    configure_file(\${path} ${GSL_BINARY_DIR}/gsl COPY_ONLY)
  endforeach ()")

add_custom_command(OUTPUT ${GSL_HEADERS}
  COMMAND ${CMAKE_COMMAND} -P copy-headers.cmake)
add_custom_target(copy-headers DEPENDS ${GSL_HEADERS})

macro(get_sources dir line source_var)
  set(${source_var})
  string(REGEX REPLACE ".*_SOURCES[ \t]*=(.*)" "\\1" sources ${line})
  string(REGEX MATCHALL "[^ \t]+" sources ${sources})
  foreach (src ${sources})
    set(${source_var} ${${source_var}} ${dir}/${src})
  endforeach ()
endmacro()

# Get subdirectories from Makefile.am.
file(STRINGS Makefile.am lines REGEX "^SUBDIRS[ \t]*=")
foreach (line ${lines})
  string(REGEX REPLACE "SUBDIRS[ \t]*=(.*)" "\\1" dirs ${line})
  string(REGEX MATCHALL "[^ ]+" dirs ${dirs})
endforeach ()

# Extract sources from automake files and add tests.
foreach (dir ${dirs})
  file(STRINGS ${dir}/Makefile.am lines)
  foreach (line ${lines})
    if (line MATCHES "_la_SOURCES[ \t]*=")
      get_sources(${dir} "${line}" SOURCES)
      if (dir STREQUAL cblas)
        add_library(gslcblas ${SOURCES})
        add_dependencies(gslcblas copy-headers)
      else ()
        set(GSL_SOURCES ${GSL_SOURCES} ${SOURCES})
      endif ()
    elseif (line MATCHES "^test.*_SOURCES[ \t]*=")
      get_sources(${dir} "${line}" SOURCES)
      string(REGEX REPLACE "(.*)_SOURCES.*" "\\1" suffix ${line})
      add_gsl_test("${dir}_${suffix}" ${SOURCES})
    endif ()
  endforeach ()
endforeach ()

add_library(gsl ${GSL_SOURCES})
add_dependencies(gsl copy-headers)
