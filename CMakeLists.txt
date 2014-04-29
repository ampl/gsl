# CMake build script for the GNU Scientific Library.

cmake_minimum_required(VERSION 2.6)
project(GSL)

set(PACKAGE "gsl")
set(PACKAGE_NAME ${PACKAGE})
set(PACKAGE_TARNAME ${PACKAGE})
set(PACKAGE_BUGREPORT "")
set(PACKAGE_URL "")
set(LT_OBJDIR ".libs/")

# The following part of config.h is hard to derive from configure.ac.
string(REGEX MATCHALL "[^\n]*\n" CONFIG
"/* Define if you have inline with C99 behavior */
#undef HAVE_C99_INLINE

/* Define to 1 if you have the declaration of `acosh', and to 0 if you don't.
   */
#undef HAVE_DECL_ACOSH

/* Define to 1 if you have the declaration of `asinh', and to 0 if you don't.
   */
#undef HAVE_DECL_ASINH

/* Define to 1 if you have the declaration of `atanh', and to 0 if you don't.
   */
#undef HAVE_DECL_ATANH

/* Define to 1 if you have the declaration of `expm1', and to 0 if you don't.
   */
#undef HAVE_DECL_EXPM1

/* Define to 1 if you have the declaration of `feenableexcept', and to 0 if
   you don't. */
#undef HAVE_DECL_FEENABLEEXCEPT

/* Define to 1 if you have the declaration of `fesettrapenable', and to 0 if
   you don't. */
#undef HAVE_DECL_FESETTRAPENABLE

/* Define to 1 if you have the declaration of `finite', and to 0 if you don't.
   */
#undef HAVE_DECL_FINITE

/* Define to 1 if you have the declaration of `fprnd_t', and to 0 if you
   don't. */
#undef HAVE_DECL_FPRND_T

/* Define to 1 if you have the declaration of `frexp', and to 0 if you don't.
   */
#undef HAVE_DECL_FREXP

/* Define to 1 if you have the declaration of `hypot', and to 0 if you don't.
   */
#undef HAVE_DECL_HYPOT

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#undef HAVE_DECL_ISFINITE

/* Define to 1 if you have the declaration of `isinf', and to 0 if you don't.
   */
#undef HAVE_DECL_ISINF

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
   */
#undef HAVE_DECL_ISNAN

/* Define to 1 if you have the declaration of `ldexp', and to 0 if you don't.
   */
#undef HAVE_DECL_LDEXP

/* Define to 1 if you have the declaration of `log1p', and to 0 if you don't.
   */
#undef HAVE_DECL_LOG1P

/* Define to 1 if you have the <dlfcn.h> header file. */
#undef HAVE_DLFCN_H

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
#undef HAVE_DOPRNT

/* Defined if you have ansi EXIT_SUCCESS and EXIT_FAILURE in stdlib.h */
#undef HAVE_EXIT_SUCCESS_AND_FAILURE

/* Defined on architectures with excess floating-point precision */
#undef HAVE_EXTENDED_PRECISION_REGISTERS

/* Define if x86 processor has sse extensions. */
#undef HAVE_FPU_X86_SSE

/* Define to 1 if you have the <ieeefp.h> header file. */
#undef HAVE_IEEEFP_H

/* Define this if IEEE comparisons work correctly (e.g. NaN != NaN) */
#undef HAVE_IEEE_COMPARISONS

/* Define this if IEEE denormalized numbers are available */
#undef HAVE_IEEE_DENORMALS

/* Define if you have inline */
#undef HAVE_INLINE

/* Define to 1 if you have the <inttypes.h> header file. */
#undef HAVE_INTTYPES_H

/* Define to 1 if you have the `m' library (-lm). */
#undef HAVE_LIBM

/* Define to 1 if you have the `memcpy' function. */
#undef HAVE_MEMCPY

/* Define to 1 if you have the `memmove' function. */
#undef HAVE_MEMMOVE

/* Define to 1 if you have the <memory.h> header file. */
#undef HAVE_MEMORY_H

/* Define this if printf can handle %Lf for long double */
#undef HAVE_PRINTF_LONGDOUBLE

/* Define to 1 if you have the <stdint.h> header file. */
#undef HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#undef HAVE_STDLIB_H

/* Define to 1 if you have the `strdup' function. */
#undef HAVE_STRDUP

/* Define to 1 if you have the <strings.h> header file. */
#undef HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#undef HAVE_STRING_H

/* Define to 1 if you have the `strtol' function. */
#undef HAVE_STRTOL

/* Define to 1 if you have the `strtoul' function. */
#undef HAVE_STRTOUL

/* Define to 1 if you have the <sys/stat.h> header file. */
#undef HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#undef HAVE_SYS_TYPES_H

/* Define to 1 if you have the <unistd.h> header file. */
#undef HAVE_UNISTD_H

/* Define to 1 if you have the `vprintf' function. */
#undef HAVE_VPRINTF

/* Define if you need to hide the static definitions of inline functions */
#undef HIDE_INLINE_STATIC

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#undef LT_OBJDIR

/* Name of package */
#undef PACKAGE

/* Define to the address where bug reports for this package should be sent. */
#undef PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#undef PACKAGE_NAME

/* Define to the full name and version of this package. */
#undef PACKAGE_STRING

/* Define to the one symbol short name of this package. */
#undef PACKAGE_TARNAME

/* Define to the home page for this package. */
#undef PACKAGE_URL

/* Define to the version of this package. */
#undef PACKAGE_VERSION

/* Defined if this is an official release */
#undef RELEASED

/* Define to 1 if you have the ANSI C header files. */
#undef STDC_HEADERS

/* Version number of package */
#undef VERSION

/* Define to 1 if type `char' is unsigned and you are not using gcc.  */
#ifndef __CHAR_UNSIGNED__
# undef __CHAR_UNSIGNED__
#endif

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#undef inline
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
#undef size_t

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
#undef volatile
")

# Get version numbers and parts of config.h from configure.ac.
file(READ configure.ac LINES)
# Replace semicolons with "<semi>" to avoid CMake messing with them.
string(REPLACE ";" "<semi>" LINES "${LINES}")
# Split into lines keeping newlines to avoid foreach skipping empty ones.
string(REGEX MATCHALL "[^\n]*\n" LINES "${LINES}")
set(ah_command FALSE)
foreach (line "${EXTRA_CONFIG}" ${LINES})
  string(REPLACE ";" "" line "${line}")
  if (ah_command)
    # Do nothing.
  elseif (line MATCHES "AC_INIT[^,]*,\\[(.*)\\].*")
    set(VERSION ${CMAKE_MATCH_1})
    message(STATUS "Got VERSION=${VERSION} from configure.ac")
  elseif (line MATCHES "((GSL|CBLAS)_(CURRENT|REVISION|AGE))=(.*)\n")
    set(${CMAKE_MATCH_1} ${CMAKE_MATCH_4})
    message(STATUS "Got ${CMAKE_MATCH_1}=${CMAKE_MATCH_4} from configure.ac")
  elseif (line MATCHES "AH_BOTTOM\\(\\[(.*)")
    set(ah_command bottom)
    set(line "${CMAKE_MATCH_1}")
  elseif (line MATCHES "AH_VERBATIM[^,]+,(.*)")
    set(ah_command verbatim)
    set(line "${CMAKE_MATCH_1}")
  endif ()
  if (ah_command)
    set(saved_ah_command ${ah_command})
    if (line MATCHES "^\\[(.*)")
      set(line "${CMAKE_MATCH_1}")
    endif ()
    if (line MATCHES "\\]\\)")
      set(ah_command FALSE)
      string(REPLACE "])" "" line "${line}")
    endif ()
    # For some reason CMake may bundle several lines together. Split them too.
    string(REGEX MATCHALL "[^\n]*\n" sublines "${line}")
    set(config_add "")
    foreach (subline ${sublines})
      set(config_add ${config_add} "${subline}")
    endforeach ()
    if (saved_ah_command STREQUAL "verbatim")
      set(CONFIG ${config_add} ${CONFIG})
    else ()
      set(CONFIG ${CONFIG} "\n" ${config_add})
    endif ()
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
#file(STRINGS ${GSL_SOURCE_DIR}/config.h.in CONFIG)
list(LENGTH CONFIG length)
math(EXPR length "${length} - 1")
foreach (i RANGE ${length})
  list(GET CONFIG ${i} line)
  if (line MATCHES "^#( *)undef (.*)\n")
    set(space "${CMAKE_MATCH_1}")
    set(var ${CMAKE_MATCH_2})
    if (NOT DEFINED ${var} OR (var MATCHES "HAVE_.*_H" AND NOT ${var}))
      set(line "/* #${space}undef ${var} */\n")
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
      set(line "#${space}define ${var} ${value}\n")
    endif ()
  endif ()
  string(REPLACE "<semi>" ";" line "${line}")
  set(CONFIG_OUT "${CONFIG_OUT}${line}")
endforeach ()
file(WRITE ${GSL_BINARY_DIR}/config.h
"/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

${CONFIG_OUT}")

include_directories(${GSL_BINARY_DIR} ${GSL_SOURCE_DIR})
add_definitions(-DHAVE_CONFIG_H)

include(CheckCCompilerFlag)
if (GSL_DISABLE_WARNINGS)
  # Disable additional warnings.
  if (MSVC)
    add_definitions(
      -D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4056 /wd4244 /wd4700 /wd4756)
  else ()
    foreach (flag -Wall -Wextra -pedantic)
      string(REPLACE ${flag} "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    endforeach ()
  endif ()
endif ()

enable_testing()

# Adds a GSL test. Usage:
#   add_gsl_test(<exename> <source> ...)
function(add_gsl_test exename)
  add_executable(${exename} ${ARGN})
  target_link_libraries(${exename} gsl gslcblas ${CMAKE_REQUIRED_LIBRARIES})
  add_test(${exename} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${exename})
endfunction()

file(GLOB GSL_HEADER_PATHS "${GSL_SOURCE_DIR}/gsl*.h")
file(GLOB dirs "${GSL_SOURCE_DIR}/*")
foreach (dir ${dirs})
  if (NOT (dir MATCHES "gsl$"))
    file(GLOB headers "${dir}/gsl*.h")
    set(GSL_HEADER_PATHS ${GSL_HEADER_PATHS} ${headers})
  endif ()
endforeach()

foreach (path ${GSL_HEADER_PATHS})
  get_filename_component(name ${path} NAME)
  set(GSL_HEADERS ${GSL_HEADERS} gsl/${name})
endforeach ()

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/copy-headers.cmake "
  file(MAKE_DIRECTORY gsl)
  foreach (path ${GSL_HEADER_PATHS})
    get_filename_component(filename \${path} NAME)
    configure_file(\${path} ${GSL_BINARY_DIR}/gsl/\${filename} COPY_ONLY)
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
foreach (dir "." ${dirs})
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

install(TARGETS gsl gslcblas
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib)

install(FILES ${GSL_HEADER_PATHS} DESTINATION include/gsl)
