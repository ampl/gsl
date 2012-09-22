# CMake build script for the GNU Scientific Library.

cmake_minimum_required(VERSION 2.6)
project(GSL)

include(CheckSymbolExists)
include(CheckCXXSourceCompiles)

# Remember to put a definition in config.h.in for each of these.

check_cxx_source_compiles("
#define _GNU_SOURCE 1
#include <fenv.h>
int main() { feenableexcept(0); }" HAVE_DECL_FEENABLEEXCEPT)

check_cxx_source_compiles("
#define _GNU_SOURCE 1
#include <fenv.h>
int main() { fesettrapenable(0); }" HAVE_DECL_FESETTRAPENABLE)

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

file(READ ${GSL_SOURCE_DIR}/config.h.in CONFIG)
foreach (var HAVE_DECL_FEENABLEEXCEPT HAVE_DECL_FESETTRAPENABLE
             HAVE_DECL_HYPOT HAVE_DECL_EXPM1 HAVE_DECL_ACOSH
             HAVE_DECL_ASINH HAVE_DECL_ATANH HAVE_DECL_LDEXP
             HAVE_DECL_FREXP HAVE_DECL_FPRND_T HAVE_DECL_ISINF
             HAVE_DECL_ISFINITE HAVE_DECL_FINITE HAVE_DECL_ISNAN)
  set(value 0)
  if (${var})
    set(value 1)
  endif ()
  string(REPLACE "#undef ${var}"
    "#define ${var} ${value}" CONFIG "${CONFIG}")
endforeach ()
file(WRITE ${GSL_BINARY_DIR}/config.h "${CONFIG}
/* Define a rounding function which moves extended precision values
   out of registers and rounds them to double-precision. This should
   be used *sparingly*, in places where it is necessary to keep
   double-precision rounding for critical expressions while running in
   extended precision. For example, the following code should ensure
   exact equality, even when extended precision registers are in use,

      double q = GSL_COERCE_DBL(3.0/7.0) ;
      if (q == GSL_COERCE_DBL(3.0/7.0)) { ... } ;

   It carries a penalty even when the program is running in double
   precision mode unless you compile a separate version of the
   library with HAVE_EXTENDED_PRECISION_REGISTERS turned off. */

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

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

#ifdef __GNUC__
#define DISCARD_POINTER(p) do { ; } while(p ? 0 : 0);
#else
#define DISCARD_POINTER(p) /* ignoring discarded pointer */
#endif

#if defined(GSL_RANGE_CHECK_OFF) || !defined(GSL_RANGE_CHECK)
#define GSL_RANGE_CHECK 0  /* turn off range checking by default internally */
#endif

#define RETURN_IF_NULL(x) if (!x) { return ; }

/* Disable deprecated functions and enums while building */
#define GSL_DISABLE_DEPRECATED 1
")

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
