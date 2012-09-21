# CMake build script for the GNU Scientific Library.

cmake_minimum_required(VERSION 2.6)
project(GSL)

include(CheckSymbolExists)
check_symbol_exists(finite "math.h" HAVE_DECL_FINITE)
check_symbol_exists(isnan "math.h" HAVE_DECL_ISNAN)

file(READ ${GSL_SOURCE_DIR}/config.h.in CONFIG)
foreach (var HAVE_DECL_FINITE HAVE_DECL_ISNAN)
  string(REPLACE "#undef ${var}"
    "#define ${var} ${${var}}" CONFIG "${CONFIG}")
endforeach ()
file(WRITE ${GSL_BINARY_DIR}/config.h "${CONFIG}
/* Substitute gsl functions for missing system functions */

#if !HAVE_DECL_HYPOT
#define hypot gsl_hypot
#endif

#if !HAVE_DECL_LOG1P
#define log1p gsl_log1p
#endif

#if !HAVE_DECL_EXPM1
#define expm1 gsl_expm1
#endif

#if !HAVE_DECL_ACOSH
#define acosh gsl_acosh
#endif

#if !HAVE_DECL_ASINH
#define asinh gsl_asinh
#endif

#if !HAVE_DECL_ATANH
#define atanh gsl_atanh
#endif

#if !HAVE_DECL_LDEXP
#define ldexp gsl_ldexp
#endif

#if !HAVE_DECL_FREXP
#define frexp gsl_frexp
#endif

#if !HAVE_DECL_ISINF
#define isinf gsl_isinf
#endif

#if !HAVE_DECL_ISFINITE
#define isfinite gsl_finite
#endif

#if !HAVE_DECL_FINITE
#define finite gsl_finite
#endif

#if !HAVE_DECL_ISNAN
#define isnan gsl_isnan
#endif
")

include_directories(${GSL_BINARY_DIR} ${GSL_SOURCE_DIR})

add_subdirectory(gsl)
add_subdirectory(utils)
add_subdirectory(sys)
add_subdirectory(test)
add_subdirectory(err)
add_subdirectory(ieee-utils)
