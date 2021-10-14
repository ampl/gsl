/*
 AMPL bindings for GNU Scientific Library.

 Copyright (C) 2012 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <math.h>
#include <stdarg.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_version.h>

#include "funcadd.h"

// Macros used for compatibility with GSL 1.x.
#if GSL_MAJOR_VERSION < 2
# define gsl_sf_mathieu_a_e gsl_sf_mathieu_a
# define gsl_sf_mathieu_b_e gsl_sf_mathieu_b
# define gsl_sf_mathieu_ce_e gsl_sf_mathieu_ce
# define gsl_sf_mathieu_se_e gsl_sf_mathieu_se
# define gsl_sf_mathieu_Mc_e gsl_sf_mathieu_Mc
# define gsl_sf_mathieu_Ms_e gsl_sf_mathieu_Ms
#endif

enum { MAX_ERROR_MESSAGE_SIZE = 100 };

static const char *const DERIVS_NOT_PROVIDED = "derivatives are not provided";

/* Computes (x / fabs(x)) * y. Returns 0 if y is 0. */
static double mul_by_sign(double x, double y) {
  return y != 0 ? (x / fabs(x)) * y : 0;
}

static char* allocate_string(arglist *al, size_t size) {
  return static_cast<char*>(al->AE->Tempmem(al->TMI, size));
}

static double* allocate_double(arglist* al, size_t size) {
  return static_cast<double*>(al->AE->Tempmem(al->TMI, size*sizeof(double)));
}
/* Formats the error message and stores it in al->Errmsg. */
static void format_error(
    arglist *al, const char *format, va_list args, char prefix) {
  size_t size = MAX_ERROR_MESSAGE_SIZE + (prefix ? 1 : 0);
  char *message = al->Errmsg = allocate_string(al, size);
  if (prefix)
    *message++ = prefix;
  al->AE->VsnprintF(message, MAX_ERROR_MESSAGE_SIZE, format, args);
}

/* Formats the error message and stores it in al->Errmsg. */
static void error(arglist *al, const char *format, ...) {
  va_list args;
  va_start(args, format);
  format_error(al, format, args, 0);
  va_end(args);
}

/* Formats the derivative error message and stores it in al->Errmsg. */
static void deriv_error(arglist *al, const char *format, ...) {
  va_list args;
  if (al->Errmsg)
    return;
  va_start(args, format);
  format_error(al, format, args, '\'');
  va_end(args);
}

/* Checks if the argument is within the bounds for derivative computation. */
static int check_deriv_arg(arglist *al, int arg, int min, int max) {
  if (arg < min) {
    deriv_error(al, "can't compute derivative: argument 'n' too small, n = %d",
                arg);
    return 0;
  }
  if (arg > max) {
    deriv_error(al, "can't compute derivative: argument 'n' too large, n = %d",
                arg);
    return 0;
  }
  return 1;
}

/* Reports a function evaluation error printing the specified suffix
 * after the function name. */
static void format_eval_error(arglist *al, char prefix, const char *suffix) {
  int n = 0, i = 0;
  size_t size = MAX_ERROR_MESSAGE_SIZE + (prefix ? 1 : 0);
  char *message = al->Errmsg = allocate_string(al, size);
  if (prefix)
    *message++ = prefix;
  n += al->AE->SnprintF(message, MAX_ERROR_MESSAGE_SIZE,
      "can't evaluate %s%s(", al->funcinfo, suffix);
  for (i = 0; i < al->n - 1; ++i) {
    n += al->AE->SnprintF(message + n, MAX_ERROR_MESSAGE_SIZE - n,
        "%g, ", al->ra[i]);
  }
  al->AE->SnprintF(message + n, MAX_ERROR_MESSAGE_SIZE - n,
      "%g)", al->ra[al->n - 1]);
}

/* Reports a function evaluation error. */
static void eval_error(arglist *al) {
  format_eval_error(al, 0, "");
}

static int check_args(arglist *al) {
  int i = 0;
  for (; i < al->n; ++i) {
    if (gsl_isnan(al->ra[i])) {
      eval_error(al);
      return 0;
    }
  }
  return 1;
}

static double check_result(arglist *al, double result) {
  int i = 0, n = 0;
  /* Function evaluation errors override derivative errors,
     because the latter may be ignored. */
  if (gsl_isnan(result)) {
    eval_error(al);
    return 0;
  }
  for (i = 0; i < al->n; ++i) {
    if (gsl_isnan(al->ra[i])) {
      eval_error(al);
      return 0;
    }
  }
  /* Check derivative errors only if there are no other errors. */
  if (al->derivs && !al->Errmsg) {
    for (i = 0; i < al->n; ++i) {
      if (gsl_isnan(al->derivs[i])) {
        format_eval_error(al, '\'', "'");
        return 0;
      }
    }
    if (al->hes) {
      for (i = 0, n = al->n * (al->n + 1) / 2; i < n; ++i) {
        if (gsl_isnan(al->hes[i])) {
          format_eval_error(al, '"', "''");
          return 0;
        }
      }
    }
  }
  return result;
}

/* Flags for check_bessel_args */
enum {
  DERIV_INT_MIN = 1 /* Derivative can be computed for n = INT_MIN */
};

/*
 * Checks whether the first argument is constant and reports an error if not.
 * Returns 1 iff the first argument is constant.
 */
static int check_const_arg(arglist *al, unsigned index, const char *name) {
  if (al->dig && al->dig[index])
    return 1;
  /* Derivative information is requested, so the argument is not constant. */
  deriv_error(al, "argument '%s' is not constant", name);
  return 0;
}

/* Checks if the argument with the specified index is representable as int. */
static int check_int_arg(arglist *al, unsigned index, const char *name) {
  double arg = al->ra[index];
  if ((int)arg != arg) {
    error(al, "argument '%s' can't be represented as int, %s = %g",
          name, name, arg);
    return 0;
  }
  if (al->derivs)
    check_const_arg(al, index, name);
  return 1;
}

/* Checks if the argument with the specified index is representable as
   unsigned int. */
static int check_uint_arg(arglist *al, unsigned index, const char *name) {
  double arg = al->ra[index];
  if ((unsigned)arg != arg) {
    error(al, "argument '%s' can't be represented as unsigned int, %s = %g",
        name, name, arg);
    return 0;
  }
  if (al->derivs)
    check_const_arg(al, index, name);
  return 1;
}

/*
 * Checks the arguments of a zero function such as gsl_sf_airy_Ai_scaled:
 * - argument with the specified index should be representable as unsigned int
 * - al->derivs should be null
 */
static int check_zero_func_args(arglist *al, unsigned s_index) {
  double arg = al->ra[s_index];
  if ((unsigned)arg != arg) {
    error(al, "argument 's' can't be represented as unsigned int, s = %g", arg);
    return 0;
  }
  if (al->derivs && check_const_arg(al, s_index, "s"))
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return 1;
}

/* Checks the arguments of a Bessel function. */
static int check_bessel_args(arglist *al, int flags, const char *arg_name) {
  int n = (int)al->ra[0];
  if (!check_int_arg(al, 0, arg_name))
    return 0;
  if (al->derivs) {
    int deriv_min = INT_MIN + ((flags & DERIV_INT_MIN) != 0 ? 0 : 1);
    if ((al->hes && !check_deriv_arg(al, n, INT_MIN + 2, INT_MAX - 2)) ||
        !check_deriv_arg(al, n, deriv_min, INT_MAX - 1)) {
      return 0;
    }
  }
  return 1;
}

#define CHECK_CALL(value, call) { \
    int status = 0; \
    gsl_sf_result result = {0, 0}; \
    status = call; \
    if (status != GSL_SUCCESS) { \
      eval_error(al); \
      return 0; \
    } \
    value = result.val; \
  }

#define ARGS1 al->ra[0]
#define ARGS2 ARGS1, al->ra[1]
#define ARGS3 ARGS2, al->ra[2]
#define ARGS4 ARGS3, al->ra[3]
#define ARGS2_PREC ARGS2, GSL_PREC_DOUBLE
#define ARGS3_PREC ARGS3, GSL_PREC_DOUBLE
#define ARGS4_PREC ARGS4, GSL_PREC_DOUBLE
#define RNG_ARGS1 rng, ARGS1
#define RNG_ARGS2 rng, ARGS2
#define RNG_ARGS3 rng, ARGS3

#define ARRAY_ARGS al->ra, 1, al->n
#define ARRAY_ARGS_AND_LAST al->ra, 1, al->n-1, al->ra[al->n-1]
#define ARRAY_ARGS_AND_LAST_TO_FIRST al->ra[al->n-1], al->ra, 1, al->n-1
#define ARRAY_ARGS_AND_LASTTWO al->ra, 1, al->n-2, al->ra[al->n-2], al->ra[al->n-1]

#define WRAP(func, args) \
  static double ampl##func(arglist *al) { \
    if (!check_args(al)) \
      return 0; \
    if (al->derivs) \
      deriv_error(al, DERIVS_NOT_PROVIDED); \
    return check_result(al, func(args)); \
  }

#define WRAP_CHECKED(func, args) \
  static double ampl##func(arglist *al) { \
    double value = 0; \
    if (!check_args(al)) \
      return 0; \
    if (al->derivs) \
      deriv_error(al, DERIVS_NOT_PROVIDED); \
    CHECK_CALL(value, func##_e(args, &result)); \
    return check_result(al, value); \
  }

#define UNUSED(x) (void)(x)

static const char *amplgsl_version(arglist *al) {
  UNUSED(al);
  return gsl_version;
}

static double amplgsl_log1p(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = 1 / (x + 1);
    if (al->hes)
      *al->hes = -deriv * deriv;
  }
  return check_result(al, gsl_log1p(x));
}

static double amplgsl_expm1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = exp(x);
    if (al->hes)
      *al->hes = deriv;
  }
  return check_result(al, gsl_expm1(x));
}

static double amplgsl_hypot(arglist *al) {
  double x = al->ra[0];
  double y = al->ra[1];
  double hypot = gsl_hypot(x, y);
  if (al->derivs) {
    real *derivs = al->derivs;
    derivs[0] = x / hypot;
    derivs[1] = y / hypot;
    if (al->hes) {
      real *hes = al->hes;
      hes[0] =  derivs[1] * derivs[1] / hypot;
      hes[1] = -derivs[0] * derivs[1] / hypot;
      hes[2] =  derivs[0] * derivs[0] / hypot;
    }
  }
  return check_result(al, hypot);
}

static double amplgsl_hypot3(arglist *al) {
  double x = al->ra[0];
  double y = al->ra[1];
  double z = al->ra[2];
  double hypot = gsl_hypot3(x, y, z);
  if (al->derivs) {
    real *derivs = al->derivs;
    derivs[0] = x / hypot;
    derivs[1] = y / hypot;
    derivs[2] = z / hypot;
    if (al->hes) {
      real *hes = al->hes;
      double dx2 = derivs[0] * derivs[0];
      double dy2 = derivs[1] * derivs[1];
      double dz2 = derivs[2] * derivs[2];
      hes[0] =  (dy2 + dz2) / hypot;
      hes[1] = -derivs[0] * derivs[1] / hypot;
      hes[2] = -derivs[0] * derivs[2] / hypot;
      hes[3] =  (dx2 + dz2) / hypot;
      hes[4] = -derivs[1] * derivs[2] / hypot;
      hes[5] =  (dx2 + dy2) / hypot;
    }
  }
  return check_result(al, hypot);
}

static double amplgsl_sf_airy_Ai(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    *al->derivs = gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
    if (al->hes)
      *al->hes = x * value;
  }
  return check_result(al, value);
}

static double amplgsl_sf_airy_Bi(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Bi(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    *al->derivs = gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
    if (al->hes)
      *al->hes = x * value;
  }
  return check_result(al, value);
}

static double amplgsl_sf_airy_Ai_scaled(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Ai_scaled(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    if (x > 0) {
      double sqrtx = sqrt(x);
      *al->derivs = gsl_sf_airy_Ai_deriv_scaled(x, GSL_PREC_DOUBLE) +
          sqrtx * gsl_sf_airy_Ai_scaled(x, GSL_PREC_DOUBLE);
      if (al->hes)
        *al->hes = (value + 4 * x * *al->derivs) / (2 * sqrtx);
    } else {
      *al->derivs = gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
      if (al->hes) {
        /* Return NaN for x = 0 since the right derivative is infinity
           and the left derivative is 0 at this point. */
        *al->hes = x != 0 ? x * value : GSL_NAN;
      }
    }
  }
  return check_result(al, value);
}

static double amplgsl_sf_airy_Bi_scaled(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Bi_scaled(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    if (x > 0) {
      double sqrtx = sqrt(x);
      *al->derivs = gsl_sf_airy_Bi_deriv_scaled(x, GSL_PREC_DOUBLE) -
          sqrtx * gsl_sf_airy_Bi_scaled(x, GSL_PREC_DOUBLE);
      if (al->hes)
        *al->hes = -(value + 4 * x * *al->derivs) / (2 * sqrtx);
    } else {
      *al->derivs = gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
      if (al->hes) {
        /* Return NaN for x = 0 since the right derivative is -infinity
           and left derivative is 0 at this point. */
        *al->hes = x != 0 ? x * value : GSL_NAN;
      }
    }
  }
  return check_result(al, value);
}

static double amplgsl_sf_airy_zero_Ai(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Ai((unsigned)al->ra[0]));
}

static double amplgsl_sf_airy_zero_Bi(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Bi((unsigned)al->ra[0]));
}

static double amplgsl_sf_airy_zero_Ai_deriv(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Ai_deriv((unsigned)al->ra[0]));
}

static double amplgsl_sf_airy_zero_Bi_deriv(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Bi_deriv((unsigned)al->ra[0]));
}

static double amplgsl_sf_bessel_J0(arglist *al) {
  double x = al->ra[0];
  double j0 = gsl_sf_bessel_J0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_J1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Jn(2, x) - j0);
  }
  return check_result(al, j0);
}

static double amplgsl_sf_bessel_J1(arglist *al) {
  double x = al->ra[0];
  double j1 = gsl_sf_bessel_J1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_J0(x) - gsl_sf_bessel_Jn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Jn(3, x) - 3 * j1);
  }
  return check_result(al, j1);
}

static double amplgsl_sf_bessel_Jn(arglist *al) {
  int n = (int)al->ra[0];
  double x = al->ra[1];
  double jn = 0;
  if (!check_bessel_args(al, 0, "n"))
    return 0;
  jn = gsl_sf_bessel_Jn(n, x);
  if (al->derivs) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Jn(n - 1, x) - gsl_sf_bessel_Jn(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Jn(n - 2, x) - 2 * jn + gsl_sf_bessel_Jn(n + 2, x));
    }
  }
  return check_result(al, jn);
}

static double amplgsl_sf_bessel_Y0(arglist *al) {
  double x = al->ra[0];
  double y0 = gsl_sf_bessel_Y0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_Y1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Yn(2, x) - y0);
  }
  return check_result(al, y0);
}

static double amplgsl_sf_bessel_Y1(arglist *al) {
  double x = al->ra[0];
  double y1 = gsl_sf_bessel_Y1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_Y0(x) - gsl_sf_bessel_Yn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Yn(3, x) - 3 * y1);
  }
  return check_result(al, y1);
}

static double amplgsl_sf_bessel_Yn(arglist *al) {
  int n = (int)al->ra[0];
  double x = al->ra[1];
  double yn = 0;
  if (!check_bessel_args(al, 0, "n"))
    return 0;
  CHECK_CALL(yn, gsl_sf_bessel_Yn_e(n, x, &result));
  if (al->derivs) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Yn(n - 1, x) - gsl_sf_bessel_Yn(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Yn(n - 2, x) - 2 * yn + gsl_sf_bessel_Yn(n + 2, x));
    }
  }
  return check_result(al, yn);
}

static double amplgsl_sf_bessel_I0(arglist *al) {
  double x = al->ra[0];
  double i0 = gsl_sf_bessel_I0(x);
  if (al->derivs) {
    *al->derivs = gsl_sf_bessel_I1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_In(2, x) + i0);
  }
  return check_result(al, i0);
}

static double amplgsl_sf_bessel_I1(arglist *al) {
  double x = al->ra[0];
  double i1 = gsl_sf_bessel_I1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_I0(x) + gsl_sf_bessel_In(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_In(3, x) + 3 * i1);
  }
  return check_result(al, i1);
}

static double amplgsl_sf_bessel_In(arglist *al) {
  int n = (int)al->ra[0];
  double x = al->ra[1];
  double in = 0;
  if (!check_bessel_args(al, 0, "n"))
    return 0;
  in = gsl_sf_bessel_In(n, x);
  if (al->derivs) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_In(n - 1, x) + gsl_sf_bessel_In(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_In(n - 2, x) + 2 * in + gsl_sf_bessel_In(n + 2, x));
    }
  }
  return check_result(al, in);
}

static double amplgsl_sf_bessel_I0_scaled(arglist *al) {
  double x = al->ra[0];
  double i0 = gsl_sf_bessel_I0_scaled(x);
  if (al->derivs) {
    double i1 = gsl_sf_bessel_I1_scaled(x);
    *al->derivs = i1 - mul_by_sign(x, i0);
    if (al->hes) {
      *al->hes = 1.5 * i0 - 2 * fabs(x) * i1 / x +
          0.5 * gsl_sf_bessel_In_scaled(2, x);
    }
  }
  return check_result(al, i0);
}

static double amplgsl_sf_bessel_I1_scaled(arglist *al) {
  double x = al->ra[0];
  double i1 = gsl_sf_bessel_I1_scaled(x);
  if (al->derivs) {
    double i0 = gsl_sf_bessel_I0_scaled(x), i2 = gsl_sf_bessel_In_scaled(2, x);
    *al->derivs = 0.5 * i0 - mul_by_sign(x, i1) + 0.5 * i2;
    if (al->hes) {
      *al->hes = -fabs(x) * i0 / x + 1.75 * i1 - fabs(x) * i2 / x +
          0.25 * gsl_sf_bessel_In_scaled(3, x);
    }
  }
  return check_result(al, i1);
}

static double amplgsl_sf_bessel_In_scaled(arglist *al) {
  int n = (int)al->ra[0];
  double x = al->ra[1];
  double in = 0;
  if (!check_bessel_args(al, 0, "n"))
    return 0;
  in = gsl_sf_bessel_In_scaled(n, x);
  if (al->derivs) {
    double in_minus_1 = gsl_sf_bessel_In_scaled(n - 1, x);
    double in_plus_1 = gsl_sf_bessel_In_scaled(n + 1, x);
    al->derivs[1] = 0.5 * in_minus_1 - mul_by_sign(x, in) + 0.5 * in_plus_1;
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_In_scaled(n - 2, x) + 6 * in +
           gsl_sf_bessel_In_scaled(n + 2, x)) -
           mul_by_sign(x, in_minus_1 + in_plus_1);
    }
  }
  return check_result(al, in);
}

static double amplgsl_sf_bessel_K0(arglist *al) {
  double x = al->ra[0];
  double k0 = gsl_sf_bessel_K0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_K1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Kn(2, x) + k0);
  }
  return check_result(al, k0);
}

static double amplgsl_sf_bessel_K1(arglist *al) {
  double x = al->ra[0];
  double k1 = gsl_sf_bessel_K1(x);
  if (al->derivs) {
    *al->derivs = -0.5 * (gsl_sf_bessel_K0(x) + gsl_sf_bessel_Kn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Kn(3, x) + 3 * k1);
  }
  return check_result(al, k1);
}

static double amplgsl_sf_bessel_Kn(arglist *al) {
  int n = (int)al->ra[0];
  double x = al->ra[1];
  double kn = 0;
  if (!check_bessel_args(al, 0, "n"))
    return 0;
  CHECK_CALL(kn, gsl_sf_bessel_Kn_e(n, x, &result));
  if (al->derivs) {
    al->derivs[1] = -0.5 *
        (gsl_sf_bessel_Kn(n - 1, x) + gsl_sf_bessel_Kn(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Kn(n - 2, x) + 2 * kn + gsl_sf_bessel_Kn(n + 2, x));
    }
  }
  return check_result(al, kn);
}

static double amplgsl_sf_bessel_K0_scaled(arglist *al) {
  double x = al->ra[0];
  double k0 = gsl_sf_bessel_K0_scaled(x);
  if (al->derivs) {
    double k1 = gsl_sf_bessel_K1_scaled(x);
    *al->derivs = k0 - k1;
    if (al->hes)
      *al->hes = 1.5 * k0 - 2 * k1 + 0.5 * gsl_sf_bessel_Kn_scaled(2, x);
  }
  return check_result(al, k0);
}

static double amplgsl_sf_bessel_K1_scaled(arglist *al) {
  double x = al->ra[0];
  double k1 = gsl_sf_bessel_K1_scaled(x);
  if (al->derivs) {
    double k0 = gsl_sf_bessel_K0_scaled(x), k2 = gsl_sf_bessel_Kn_scaled(2, x);
    *al->derivs = -0.5 * k0 + k1 - 0.5 * k2;
    if (al->hes)
      *al->hes = -k0 + 1.75 * k1 - k2 + 0.25 * gsl_sf_bessel_Kn_scaled(3, x);
  }
  return check_result(al, k1);
}

static double amplgsl_sf_bessel_Kn_scaled(arglist *al) {
  int n = (int)al->ra[0];
  double x = al->ra[1];
  double kn = 0;
  if (!check_bessel_args(al, 0, "n"))
    return 0;
  CHECK_CALL(kn, gsl_sf_bessel_Kn_scaled_e(n, x, &result));
  if (al->derivs) {
    double kn_minus_1 = gsl_sf_bessel_Kn_scaled(n - 1, x);
    double kn_plus_1 = gsl_sf_bessel_Kn_scaled(n + 1, x);
    al->derivs[1] = -0.5 * (kn_minus_1 - 2 * kn + kn_plus_1);
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Kn_scaled(n - 2, x) - 4 * kn_minus_1 + 6 * kn -
              4 * kn_plus_1 + gsl_sf_bessel_Kn_scaled(n + 2, x));
    }
  }
  return check_result(al, kn);
}

static double amplgsl_sf_bessel_j0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? (x * cos(x) - sin(x)) / gsl_pow_2(x) : 0;
    if (al->hes)
      *al->hes = x != 0 ?
          ((2 - gsl_pow_2(x)) * sin(x) - 2 * x * cos(x)) / gsl_pow_3(x) :
          -1.0 / 3.0;
  }
  return check_result(al, gsl_sf_bessel_j0(x));
}

static double amplgsl_sf_bessel_j1(arglist *al) {
  double x = al->ra[0];
  double j1 = gsl_sf_bessel_j1(x);
  if (al->derivs) {
    *al->derivs = x != 0 ? (sin(x) - 2 * j1) / x : 1.0 / 3.0;
    if (al->hes) {
      *al->hes = x != 0 ? (x * (gsl_pow_2(x) - 6) * cos(x) -
          3 * (gsl_pow_2(x) - 2) * sin(x)) / gsl_pow_4(x) : 0;
    }
  }
  return check_result(al, j1);
}

static double amplgsl_sf_bessel_j2(arglist *al) {
  double x = al->ra[0];
  double j2 = gsl_sf_bessel_j2(x);
  if (al->derivs) {
    *al->derivs = x != 0 ? gsl_sf_bessel_j1(x) - 3 * j2 / x : 0;
    if (al->hes) {
      *al->hes = x != 0 ? (x * (5 * gsl_pow_2(x) - 36) * cos(x) +
          (gsl_pow_4(x) - 17 * gsl_pow_2(x) + 36) * sin(x)) / gsl_pow_5(x) :
              2.0 / 15.0;
    }
  }
  return check_result(al, j2);
}

static double amplgsl_sf_bessel_jl(arglist *al) {
  int el = (int)al->ra[0];
  double x = al->ra[1];
  double jl = 0;
  if (!check_bessel_args(al, DERIV_INT_MIN, "l"))
    return 0;
  jl = gsl_sf_bessel_jl(el, x);
  if (al->derivs) {
    double jl_plus_1 = gsl_sf_bessel_jl(el + 1, x);
    if (x == 0)
      al->derivs[1] = el == 1 ? 1.0 / 3.0 : 0;
    else
      al->derivs[1] = el * jl / x - jl_plus_1;
    if (al->hes) {
      double hes = 0;
      if (x == 0) {
        if (el == 0)
          hes = -1.0 / 3.0;
        else if (el == 2)
          hes = 2.0 / 15.0;
      } else {
        double jl_minus_1 = 0, jl_minus_2 = 0;
        if (el == 0) {
          jl_minus_1 = cos(x) / x;
          jl_minus_2 = -(cos(x) / x + sin(x)) / x;
        } else if (el == 1) {
          jl_minus_1 = gsl_sf_bessel_jl(el - 1, x);
          jl_minus_2 = cos(x) / x;
        } else {
          jl_minus_1 = gsl_sf_bessel_jl(el - 1, x);
          jl_minus_2 = gsl_sf_bessel_jl(el - 2, x);
        }
        hes = (
          gsl_pow_2(x) * jl_minus_2 -
          2 * x * jl_minus_1 -
          (2 * gsl_pow_2(x) - 3) * jl + 2 * x * jl_plus_1 +
          gsl_pow_2(x) * gsl_sf_bessel_jl(el + 2, x)) / (4 * gsl_pow_2(x));
      }
      al->hes[2] = hes;
    }
  }
  return check_result(al, jl);
}

static double amplgsl_sf_bessel_y0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = (x * sin(x) + cos(x)) / gsl_pow_2(x);
    if (al->hes) {
      *al->hes = ((gsl_pow_2(x) - 2) * cos(x) - 2 * x * sin(x)) /
          gsl_pow_3(x);
    }
  }
  return check_result(al, gsl_sf_bessel_y0(x));
}

static double amplgsl_sf_bessel_y1(arglist *al) {
  double x = al->ra[0];
  double y1 = gsl_sf_bessel_y1(x);
  if (al->derivs) {
    *al->derivs = -(2 * y1 + cos(x)) / x;
    if (al->hes) {
      *al->hes = (x * (gsl_pow_2(x) - 6) * sin(x) +
          3 * (gsl_pow_2(x) - 2) * cos(x)) / gsl_pow_4(x);
    }
  }
  return check_result(al, y1);
}

static double amplgsl_sf_bessel_y2(arglist *al) {
  double x = al->ra[0];
  double y2 = gsl_sf_bessel_y2(x);
  if (al->derivs) {
    double y1 = gsl_sf_bessel_y1(x);
    *al->derivs = y1 - (3 * y2) / x;
    if (al->hes) {
      *al->hes = ((36 - 5 * gsl_pow_2(x)) * y1 -
          (gsl_pow_2(x) - 12) * cos(x)) / gsl_pow_3(x);
    }
  }
  return check_result(al, y2);
}

static double amplgsl_sf_bessel_yl(arglist *al) {
  int el = (int)al->ra[0];
  double x = al->ra[1];
  double yl = 0;
  if (!check_bessel_args(al, 0, "l"))
    return 0;
  yl = gsl_sf_bessel_yl(el, x);
  if (al->derivs) {
    double yl_minus_1 = el != 0 ? gsl_sf_bessel_yl(el - 1, x) : sin(x) / x;
    double yl_plus_1 = gsl_sf_bessel_yl(el + 1, x);
    al->derivs[1] = 0.5 * (yl_minus_1 - yl / x - yl_plus_1);
    if (al->hes) {
      double yl_minus_2 = 0;
      if (el == 0)
        yl_minus_2 = (cos(x) - sin(x) / x) / x;
      else if (el == 1)
        yl_minus_2 = sin(x) / x;
      else
        yl_minus_2 = gsl_sf_bessel_yl(el - 2, x);
      al->hes[2] = (
          gsl_pow_2(x) * yl_minus_2 - 2 * x * yl_minus_1 -
          (2 * gsl_pow_2(x) - 3) * yl + 2 * x * yl_plus_1 +
          gsl_pow_2(x) * gsl_sf_bessel_yl(el + 2, x)) / (4 * gsl_pow_2(x));
    }
  }
  return check_result(al, yl);
}

static double amplgsl_sf_bessel_i0_scaled(arglist *al) {
  double x = al->ra[0];
  double i0 = gsl_sf_bessel_i0_scaled(x);
  if (al->derivs) {
    /* Contrary to the documentation, gsl_sf_bessel_i0_scaled
       implements \exp(-|x|) \sqrt{\pi}/\sqrt{2x} I_{1/2}(x)
       and not \exp(-|x|) \sqrt{\pi/(2x)} I_{1/2}(x).
       These are different since \sqrt(1/x) != \sqrt(x) for negative x. */
    double hyp_coef = exp(-fabs(x)) / x;
    double i_minus_1 = hyp_coef * cosh(x);
    double i1 = gsl_sf_bessel_i1_scaled(x);
    double coef = -(1 + 2 * fabs(x)) / x;
    *al->derivs = 0.5 * (i_minus_1 + coef * i0 + i1);
    if (al->hes) {
      coef *= 2;
      *al->hes = 0.25 * (
          hyp_coef * sinh(x) - i_minus_1 / x +
          coef * i_minus_1 +
          (3 + 6 * gsl_pow_2(x) + 4 * fabs(x)) * i0 / gsl_pow_2(x) +
          coef * i1 +
          gsl_sf_bessel_il_scaled(2, x));
    }
  }
  return check_result(al, i0);
}

static double amplgsl_sf_bessel_i1_scaled(arglist *al) {
  double x = al->ra[0];
  double i1 = gsl_sf_bessel_i1_scaled(x);
  if (al->derivs) {
    /* Contrary to the documentation, gsl_sf_bessel_i1_scaled
       implements \exp(-|x|) \sqrt{\pi}/\sqrt{2x} I_{1+1/2}(x)
       and not \exp(-|x|) \sqrt{\pi/(2x)} I_{1+1/2}(x).
       These are different since \sqrt(1/x) != \sqrt(x) for negative x. */
    double i0 = gsl_sf_bessel_i0_scaled(x);
    double i2 = gsl_sf_bessel_i2_scaled(x);
    double coef = -(1 + 2 * fabs(x)) / x;
    *al->derivs = x != 0 ? 0.5 * (i0 + coef * i1 + i2) : 1.0 / 3.0;
    if (al->hes) {
      coef *= 2;
      *al->hes = 0.25 * (
          exp(-fabs(x)) * cosh(x) / x +
          coef * i0 +
          (3 + 6 * gsl_pow_2(x) + 4 * fabs(x)) * i1 / gsl_pow_2(x) +
          coef * i2 +
          gsl_sf_bessel_il_scaled(3, x));
    }
  }
  return check_result(al, i1);
}

static double amplgsl_sf_bessel_i2_scaled(arglist *al) {
  double x = al->ra[0];
  double i2 = gsl_sf_bessel_i2_scaled(x);
  if (al->derivs) {
    /* Contrary to the documentation, gsl_sf_bessel_i2_scaled
       implements \exp(-|x|) \sqrt{\pi}/\sqrt{2x} I_{2+1/2}(x)
       and not \exp(-|x|) \sqrt{\pi/(2x)} I_{2+1/2}(x).
       These are different since \sqrt(1/x) != \sqrt(x) for negative x. */
    double i1 = gsl_sf_bessel_i1_scaled(x);
    double i3 = gsl_sf_bessel_il_scaled(3, x);
    double coef = -(1 + 2 * fabs(x)) / x;
    *al->derivs = x != 0 ? 0.5 * (i1 + coef * i2 + i3) : 0;
    if (al->hes) {
      coef *= 2;
      *al->hes = x != 0 ? 0.25 * (
          gsl_sf_bessel_i0_scaled(x) +
          coef * i1 +
          (3 + 6 * gsl_pow_2(x) + 4 * fabs(x)) * i2 / gsl_pow_2(x) +
          coef * i3 +
          gsl_sf_bessel_il_scaled(4, x)) : 2.0 / 15.0;
    }
  }
  return check_result(al, i2);
}

static double amplgsl_sf_bessel_il_scaled(arglist *al) {
  int el = (int)al->ra[0];
  double x = al->ra[1];
  double il = 0;
  if (!check_bessel_args(al, 0, "l"))
    return 0;
  il = gsl_sf_bessel_il_scaled(el, x);
  if (al->derivs) {
    double il_minus_1 = el != 0 ?
        gsl_sf_bessel_il_scaled(el - 1, x) : exp(-fabs(x)) * cosh(x) / x;
    double il_plus_1 = gsl_sf_bessel_il_scaled(el + 1, x);
    double coef = -(1 + 2 * fabs(x)) / x;
    double deriv = GSL_NAN;
    if (x == 0) {
      /* If el <= 0, keep deriv equal to GSL_NAN. */
      if (el == 1)
        deriv = 1.0 / 3.0;
      else if (el > 1)
        deriv = 0;
    } else {
      deriv = 0.5 * (il_minus_1 + coef * il + il_plus_1);
    }
    al->derivs[1] = deriv;
    if (al->hes) {
      double hes = GSL_NAN;
      double il_minus_2 = 0;
      if (el == 0)
        il_minus_2 = (exp(-fabs(x)) * (sinh(x) - cosh(x) / x)) / x;
      else if (el == 1)
        il_minus_2 = exp(-fabs(x)) * cosh(x) / x;
      else
        il_minus_2 = gsl_sf_bessel_il_scaled(el - 2, x);
      coef *= 2;
      if (x == 0) {
        /* If el == 1 or el < 0, keep hes equal to GSL_NAN. */
        if (el == 0)
          hes = 4.0 / 3.0;
        else if (el == 2)
          hes = 2.0 / 15.0;
        else if (el > 2)
          hes = 0;
      } else {
        hes = 0.25 * (
          il_minus_2 +
          coef * il_minus_1 +
          (3 + 4 * fabs(x) + 6 * gsl_pow_2(x)) * il / gsl_pow_2(x) +
          coef * il_plus_1 +
          gsl_sf_bessel_il_scaled(el + 2, x));
      }
      al->hes[2] = hes;
    }
  }
  return check_result(al, il);
}

static double amplgsl_sf_bessel_k0_scaled(arglist *al) {
  double x = al->ra[0];
  double k0 = gsl_sf_bessel_k0_scaled(x);
  if (al->derivs) {
    double pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -pi_sqrt_inv_x / (2 * pow(x, 1.5));
    if (al->hes)
      *al->hes = pi_sqrt_inv_x / pow(x, 2.5);
  }
  return check_result(al, k0);
}

static double amplgsl_sf_bessel_k1_scaled(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -(pi_sqrt_inv_x * (x + 2)) / (2 * pow(x, 2.5));
    if (al->hes)
      *al->hes = (pi_sqrt_inv_x * (x + 3)) / pow(x, 3.5);
  }
  return check_result(al, gsl_sf_bessel_k1_scaled(x));
}

static double amplgsl_sf_bessel_k2_scaled(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -pi_sqrt_inv_x * (x + 3) * (x + 3) / (2 * pow(x, 3.5));
    if (al->hes)
      *al->hes = pi_sqrt_inv_x * (x * x + 9 * x + 18) / pow(x, 4.5);
  }
  return check_result(al, gsl_sf_bessel_k2_scaled(x));
}

static double amplgsl_sf_bessel_kl_scaled(arglist *al) {
  int el = (int)al->ra[0];
  double x = al->ra[1];
  double kl = 0;
  if (!check_bessel_args(al, 0, "l"))
    return 0;
  kl = gsl_sf_bessel_kl_scaled(el, x);
  if (al->derivs) {
    double kl_minus_1 = el != 0 ?
        gsl_sf_bessel_kl_scaled(el - 1, x) : M_PI_2 / x;
    double kl_plus_1 = gsl_sf_bessel_kl_scaled(el + 1, x);
    double coef = (1 - 2 * x) / x;
    al->derivs[1] = -0.5 * (kl_minus_1 + coef * kl + kl_plus_1);
    if (al->hes) {
      double kl_minus_2 = 0;
      if (el == 0)
        kl_minus_2 = M_PI_2 * (1 / x + 1) / x;
      else if (el == 1)
        kl_minus_2 = M_PI_2 / x;
      else
        kl_minus_2 = gsl_sf_bessel_kl_scaled(el - 2, x);
      coef *= 2;
      al->hes[2] = 0.25 * (
          kl_minus_2 +
          coef * kl_minus_1 +
          (3 - 4 * x + 6 * gsl_pow_2(x)) * kl / gsl_pow_2(x) +
          coef * kl_plus_1 +
          gsl_sf_bessel_kl_scaled(el + 2, x));
    }
  }
  return check_result(al, kl);
}

static double amplgsl_sf_bessel_Jnu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double jn = gsl_sf_bessel_Jnu(n, x);
  if (al->derivs && check_const_arg(al, 0, "nu")) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Jnu(n - 1, x) - gsl_sf_bessel_Jnu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Jnu(n - 2, x) - 2 * jn + gsl_sf_bessel_Jnu(n + 2, x));
    }
  }
  return check_result(al, jn);
}

static double amplgsl_sf_bessel_Ynu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double yn = gsl_sf_bessel_Ynu(n, x);
  if (al->derivs && check_const_arg(al, 0, "nu")) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Ynu(n - 1, x) - gsl_sf_bessel_Ynu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Ynu(n - 2, x) - 2 * yn + gsl_sf_bessel_Ynu(n + 2, x));
    }
  }
  return check_result(al, yn);
}

static double amplgsl_sf_bessel_Inu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double in = gsl_sf_bessel_Inu(n, x);
  if (al->derivs && check_const_arg(al, 0, "nu")) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Inu(n - 1, x) + gsl_sf_bessel_Inu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Inu(n - 2, x) + 2 * in + gsl_sf_bessel_Inu(n + 2, x));
    }
  }
  return check_result(al, in);
}

static double amplgsl_sf_bessel_Inu_scaled(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double in = gsl_sf_bessel_Inu_scaled(n, x);
  if (al->derivs && check_const_arg(al, 0, "nu")) {
    double in_minus_1 = 0, in_plus_1 = 0;
    in_minus_1 = gsl_sf_bessel_Inu_scaled(n - 1, x);
    in_plus_1 = gsl_sf_bessel_Inu_scaled(n + 1, x);
    al->derivs[1] = 0.5 * in_minus_1 - mul_by_sign(x, in) + 0.5 * in_plus_1;
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Inu_scaled(n - 2, x) + 6 * in +
              gsl_sf_bessel_Inu_scaled(n + 2, x)) -
              mul_by_sign(x, in_minus_1 + in_plus_1);
    }
  }
  return check_result(al, in);
}

static double amplgsl_sf_bessel_Knu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double kn = gsl_sf_bessel_Knu(n, x);
  if (al->derivs && check_const_arg(al, 0, "nu")) {
    al->derivs[1] = -0.5 *
        (gsl_sf_bessel_Knu(n - 1, x) + gsl_sf_bessel_Knu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Knu(n - 2, x) + 2 * kn + gsl_sf_bessel_Knu(n + 2, x));
    }
  }
  return check_result(al, kn);
}

static double amplgsl_sf_bessel_lnKnu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  if (al->derivs && check_const_arg(al, 0, "nu")) {
    double kn = 0, kn_minus_1_plus_1 = 0;
    kn = gsl_sf_bessel_Knu(n, x);
    kn_minus_1_plus_1 =
        gsl_sf_bessel_Knu(n - 1, x) + gsl_sf_bessel_Knu(n + 1, x);
    al->derivs[1] = -0.5 * kn_minus_1_plus_1 / kn;
    if (al->hes) {
      al->hes[2] = 0.25 *
          (kn * (gsl_sf_bessel_Knu(n - 2, x) + 2 * kn +
          gsl_sf_bessel_Knu(n + 2, x)) -
              kn_minus_1_plus_1 * kn_minus_1_plus_1) / (kn * kn);
    }
  }
  return check_result(al, gsl_sf_bessel_lnKnu(n, x));
}

static double amplgsl_sf_bessel_Knu_scaled(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double kn = gsl_sf_bessel_Knu_scaled(n, x);
  if (al->derivs && check_const_arg(al, 0, "nu")) {
    double kn_minus_1 = 0, kn_plus_1 = 0;
    kn_minus_1 = gsl_sf_bessel_Knu_scaled(n - 1, x);
    kn_plus_1 = gsl_sf_bessel_Knu_scaled(n + 1, x);
    al->derivs[1] = -0.5 * (kn_minus_1 - 2 * kn + kn_plus_1);
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Knu_scaled(n - 2, x) - 4 * kn_minus_1 + 6 * kn -
              4 * kn_plus_1 + gsl_sf_bessel_Knu_scaled(n + 2, x));
    }
  }
  return check_result(al, kn);
}

static double amplgsl_sf_bessel_zero_J0(arglist *al) {
  double value = 0;
  if (!check_zero_func_args(al, 0))
    return 0;
  CHECK_CALL(value, gsl_sf_bessel_zero_J0_e((unsigned)al->ra[0], &result));
  return check_result(al, value);
}

static double amplgsl_sf_bessel_zero_J1(arglist *al) {
  double value = 0;
  if (!check_zero_func_args(al, 0))
    return 0;
  CHECK_CALL(value, gsl_sf_bessel_zero_J1_e((unsigned)al->ra[0], &result));
  return check_result(al, value);
}

static double amplgsl_sf_bessel_zero_Jnu(arglist *al) {
  double value = 0;
  double nu = al->ra[0];
  if (!check_zero_func_args(al, 1))
    return 0;
  CHECK_CALL(value, gsl_sf_bessel_zero_Jnu_e(nu, (unsigned)al->ra[1], &result));
  return check_result(al, value);
}

static double amplgsl_sf_clausen(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -log(2 * sin(0.5 * fabs(x)));
    if (al->hes)
      *al->hes = fmod(x, M_PI) != 0 ? -0.5 * tan(M_PI_2 - 0.5 * x) : GSL_NAN;
  }
  return check_result(al, gsl_sf_clausen(x));
}

static double amplgsl_sf_hydrogenicR_1(arglist *al) {
  double Z = al->ra[0], r = al->ra[1];
  if (al->derivs) {
    real *derivs = al->derivs;
    double exp_minusZr = exp(-Z * r);
    derivs[0] = sqrt(Z) * exp_minusZr * (3 - 2 * r * Z);
    derivs[1] = -2 * pow(Z, 2.5) * exp_minusZr;
    if (al->hes) {
      real *hes = al->hes;
      hes[0] = (exp_minusZr * (4 * gsl_pow_2(r * Z) - 12 * r * Z + 3)) /
          (2 *sqrt(Z));
      hes[1] = pow(Z, 1.5) * exp_minusZr * (2 * r * Z - 5);
      hes[2] = 2 * pow(Z, 3.5) * exp_minusZr;
    }
  }
  return check_result(al, gsl_sf_hydrogenicR_1(Z, r));
}

static double amplgsl_sf_hydrogenicR(arglist *al) {
  double value = 0;
  if (!check_int_arg(al, 0, "n") || !check_int_arg(al, 1, "l"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  CHECK_CALL(value, gsl_sf_hydrogenicR_e(
      (int)al->ra[0], (int)al->ra[1], al->ra[2], al->ra[3], &result))
  return check_result(al, value);
}

static double amplgsl_sf_coulomb_CL(arglist *al) {
  gsl_sf_result result = {0, 0};
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  if (gsl_sf_coulomb_CL_e(al->ra[0], al->ra[1], &result)) {
    eval_error(al);
    return 0;
  }
  return check_result(al, result.val);
}

static int check_coupling_args(arglist *al, const char *const* arg_names) {
  unsigned i = 0, n_args = al->n;
  for (; i < n_args; ++i) {
    if (!check_int_arg(al, i, arg_names[i]))
      return 0;
  }
  return 1;
}

static double amplgsl_sf_coupling_3j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_ma", "two_mb", "two_mc"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return check_result(al, gsl_sf_coupling_3j(
      (int)al->ra[0], (int)al->ra[1], (int)al->ra[2],
      (int)al->ra[3], (int)al->ra[4], (int)al->ra[5]));
}

static double amplgsl_sf_coupling_6j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_jd", "two_je", "two_jf"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return check_result(al, gsl_sf_coupling_6j(
      (int)al->ra[0], (int)al->ra[1], (int)al->ra[2],
      (int)al->ra[3], (int)al->ra[4], (int)al->ra[5]));
}

static double amplgsl_sf_coupling_9j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_jd", "two_je", "two_jf",
      "two_jg", "two_jh", "two_ji"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return check_result(al, gsl_sf_coupling_9j(
      (int)al->ra[0], (int)al->ra[1], (int)al->ra[2],
      (int)al->ra[3], (int)al->ra[4], (int)al->ra[5],
      (int)al->ra[6], (int)al->ra[7], (int)al->ra[8]));
}

static double amplgsl_sf_dawson(arglist *al) {
  double x = al->ra[0];
  double f = gsl_sf_dawson(x);
  if (al->derivs) {
    double deriv = *al->derivs = 1 - 2 * x * f;
    if (al->hes)
      *al->hes = - 2 * (f + x * deriv);
  }
  return check_result(al, f);
}

/* Values of the right derivatives of the Debye functions at 0. */
static const double DEBYE_DERIV_AT_0[] = {
    -1.0 / 4.0, -1.0 / 3.0,  -3.0 / 8.0,
    -2.0 / 5.0, -5.0 / 12.0, -3.0 / 7.0
};

/* Values of the second right derivatives of the Debye functions at 0. */
static const double DEBYE_DERIV2_AT_0[] = {
    1.0 / 18.0, 1.0 / 12.0, 1.0 / 10.0,
    1.0 / 9.0,  5.0 / 42.0, 1.0 / 8.0
};

static double debye(arglist *al, int n, double (*func)(double)) {
  double x = al->ra[0];
  double f = func(x);
  if (al->derivs) {
    double exp_x = exp(x);
    double deriv = *al->derivs = x != 0 ?
        n * (1 / (exp_x - 1) - f / x) : DEBYE_DERIV_AT_0[n - 1];
    if (al->hes) {
      *al->hes = x != 0 ? n * (-exp_x / gsl_pow_2(exp_x - 1) +
          f / gsl_pow_2(x) - deriv / x) : DEBYE_DERIV2_AT_0[n - 1];
    }
  }
  return check_result(al, f);
}

#define DEBYE(n) \
  static double amplgsl_sf_debye_##n(arglist *al) { \
    return debye(al, n, gsl_sf_debye_##n); \
  }

DEBYE(1)
DEBYE(2)
DEBYE(3)
DEBYE(4)
DEBYE(5)
DEBYE(6)

static double amplgsl_sf_dilog(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = 0;
    if (x == 0) {
      deriv = 1;
    } else if (x == 1) {
      deriv = GSL_POSINF;
    } else {
      gsl_complex log = gsl_complex_log(gsl_complex_rect(1 - x, 0));
      deriv = -GSL_REAL(log) / x;
    }
    *al->derivs = deriv;
    if (al->hes)
      *al->hes = x != 0 ? (1 / (1 - x) - deriv) / x : 0.5;
  }
  return check_result(al, gsl_sf_dilog(x));
}

static double amplgsl_sf_ellint_Kcomp(arglist *al) {
  double k = al->ra[0];
  double kcomp = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double ecomp = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
    double divisor = k * (1 - k * k);
    *al->derivs = k != 0 ? ecomp / divisor - kcomp / k : 0;
    if (al->hes) {
      *al->hes = k != 0 ?
          ((2 * gsl_pow_4(k) - 3 * gsl_pow_2(k) + 1) * kcomp +
          (3 * gsl_pow_2(k) - 1) * ecomp) / gsl_pow_2(divisor) : M_PI_4;
    }
  }
  return check_result(al, kcomp);
}

static double amplgsl_sf_ellint_Ecomp(arglist *al) {
  double k = al->ra[0];
  double ecomp = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double kcomp = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
    *al->derivs = k != 0 ? (ecomp - kcomp) / k : 0;
    if (al->hes) {
      *al->hes = k != 0 ?
          ((k * k - 1) * kcomp + ecomp) / (k * k * (k * k - 1)) : -M_PI_4;
    }
  }
  return check_result(al, ecomp);
}

static double amplgsl_sf_ellint_Pcomp(arglist *al) {
  double k = al->ra[0], n = al->ra[1];
  double pcomp = gsl_sf_ellint_Pcomp(k, n, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double ecomp = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
    double kcomp = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
    double divisor = (k * k - 1) * (k * k + n);
    if (k != 0 || n != 0) {
      al->derivs[0] = -k * ((k * k - 1) * pcomp + ecomp) / divisor;
      if (n != 0) {
        al->derivs[1] = (-kcomp * (k * k + n) +
          (k * k - n * n) * pcomp + n * ecomp) /
          (2 * n * (n + 1) * (k * k + n));
      } else {
        al->derivs[1] =
            -(4 * kcomp + M_PI * k * k * gsl_sf_hyperg_2F1(0.5, 1.5, 2, k * k) -
                4 * ecomp) / (8 * k * k);
      }
    } else {
      al->derivs[0] = 0;
      al->derivs[1] = -M_PI_4;
    }
    if (al->hes) {
      if (k != 0 || n != 0) {
        al->hes[0] = ((k * k - 1) * (kcomp * (k * k + n) +
            (k * k - 1) * (2 * k * k - n) * pcomp) +
            (3 * gsl_pow_4(k) - k * k + 2 * n) * ecomp) / gsl_pow_2(divisor);
        al->hes[1] = (k * ((k * k - 1) * (kcomp * (k * k + n) +
            (n * (3 * n + 2) - k * k) * pcomp) +
            n * (-k * k + 2 * n + 3) * ecomp)) /
            (2 * divisor * n * (n + 1) * (k * k + n));
        al->hes[2] = (kcomp * (gsl_pow_4(k) * (4 * n + 1) +
            3 * k * k * n * (3 * n + 1) + n * n * (5 * n + 2)) +
            n * (k * k * (1 - 2 * n) - n * (5 * n + 2)) * ecomp -
            (gsl_pow_4(k) * (4 * n + 1) + 2 * k * k * n * (5 * n + 2) -
                3 * gsl_pow_4(n)) * pcomp) /
                (4 * gsl_pow_2(n * (n + 1) * (k * k + n)));
      } else {
        al->hes[0] = M_PI_4;
        al->hes[1] = 0;
        al->hes[2] = 3 * M_PI / 8;
      }
    }
  }
  return check_result(al, pcomp);
}

static double amplgsl_sf_ellint_F(arglist *al) {
  double phi = al->ra[0], k = al->ra[1];
  double f = gsl_sf_ellint_F(phi, k, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double e = gsl_sf_ellint_E(phi, k, GSL_PREC_DOUBLE);
    al->derivs[0] = 1 / sqrt(1 - gsl_pow_2(k * sin(phi)));
    if (!al->dig || !al->dig[1]) {
      if (k == 0) {
        al->derivs[1] = 0;
      } else if (fabs(k) == 1) {
        double sec_phi = 1 / cos(phi);
        al->derivs[1] = 0.5 * k * (gsl_atanh(sin(phi)) -
           sec_phi * sec_phi * ((1 + cos(2 * phi)) * log(sec_phi + tan(phi)) -
               sin(phi)));
      } else {
        al->derivs[1] = (e + (k * k - 1) * f -
            (k * k * cos(phi) * sin(phi)) / sqrt(1 - gsl_pow_2(k * sin(phi)))) /
            (k - gsl_pow_3(k));
      }
    }
    if (al->hes) {
      double k2 = k * k;
      al->hes[0] = (k2 * sin(phi) * cos(phi)) /
          pow(1 - gsl_pow_2(k * sin(phi)), 1.5);
      al->hes[1] = k * gsl_pow_2(sin(phi)) /
          pow(1 - gsl_pow_2(k * sin(phi)), 1.5);
      if (k == 0) {
        al->hes[2] = 0.5 * (phi - cos(phi) * sin(phi));
      } else if (fabs(k) == 1) {
        double sec_phi = 1 / cos(phi);
        al->hes[2] = sec_phi * sec_phi *
            (gsl_atanh(sin(phi)) * (2 - 46 * cos(2 * phi)) +
                8 * (1 + 7 * cos(2 * phi)) * log(sec_phi + tan(phi)) +
                sec_phi * (-11 * sec_phi * sin(3 * phi) + 13 * tan(phi))) / 32;
      } else {
        /* sub1 and sub2 are just common subexpressions */
        double sub1 = 1 - 3 * k2;
        double sub2 = M_SQRT2 * pow(2 - k2 + k2 * cos(2 * phi), 1.5);
        al->hes[2] = -(sub1 * sub2 * e - (sub1 + 2 * gsl_pow_4(k)) * sub2 * f +
          4 * gsl_pow_4(k) * (sub1 * cos(phi) * gsl_pow_3(sin(phi)) +
              sin(2 * phi))) / (gsl_pow_2(k * (k2 - 1)) * sub2);
      }
    }
  }
  return check_result(al, f);
}

static double amplgsl_sf_ellint_E(arglist *al) {
  double phi = al->ra[0], k = al->ra[1];
  double e = gsl_sf_ellint_E(phi, k, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double f = gsl_sf_ellint_F(phi, k, GSL_PREC_DOUBLE);
    double d_phi = al->derivs[0] = sqrt(1 - gsl_pow_2(k * sin(phi)));
    al->derivs[1] = k != 0 ? (e - f) / k : 0;
    if (al->hes) {
      double k2 = k * k;
      al->hes[0] = -k2 * cos(phi) * sin(phi) / d_phi;
      al->hes[1] = -k * gsl_pow_2(sin(phi)) / d_phi;
      if (k == 0) {
        al->hes[2] = -0.5 * phi + 0.25 * sin(2 * phi);
      } else if (fabs(k) == 1) {
        double sec_phi = 1 / cos(phi), tan_phi = tan(phi);
        al->hes[2] = -0.5 * gsl_atanh(sin(phi)) + log(sec_phi + tan_phi) -
            0.5 * sec_phi * tan_phi;
      } else {
        al->hes[2] = ((k2 - 1) * sqrt(4 - 2 * k2 + 2 * k2 * cos(2 * phi)) * f +
          2 * e * d_phi - k2 * sin(2 * phi)) / (2 * k2 * (k2 - 1) * d_phi);
      }
    }
  }
  return check_result(al, e);
}

WRAP_CHECKED(gsl_sf_ellint_P, ARGS3_PREC)

#if GSL_MAJOR_VERSION >= 2
# define GSL_ELLINT_D_ARGS ARGS2_PREC
#else
# define GSL_ELLINT_D_ARGS ARGS2, 0, GSL_PREC_DOUBLE
#endif
WRAP_CHECKED(gsl_sf_ellint_D, GSL_ELLINT_D_ARGS)

WRAP_CHECKED(gsl_sf_ellint_RC, ARGS2_PREC)
WRAP_CHECKED(gsl_sf_ellint_RD, ARGS3_PREC)
WRAP_CHECKED(gsl_sf_ellint_RF, ARGS3_PREC)
WRAP_CHECKED(gsl_sf_ellint_RJ, ARGS4_PREC)

static double amplgsl_sf_erf(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 2 * exp(-x * x) / sqrt(M_PI);
    if (al->hes)
      *al->hes = -2 * x * *al->derivs;
  }
  return check_result(al, gsl_sf_erf(x));
}

static double amplgsl_sf_erfc(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -2 * exp(-x * x) / sqrt(M_PI);
    if (al->hes)
      *al->hes = -2 * x * *al->derivs;
  }
  return check_result(al, gsl_sf_erfc(x));
}

static double amplgsl_sf_log_erfc(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double erfc = gsl_sf_erfc(x);
    *al->derivs = -2 * exp(-x * x) / (sqrt(M_PI) * erfc);
    if (al->hes) {
      *al->hes = -2 * x * *al->derivs -
          ((4 * exp(-2 * x * x)) / (M_PI * erfc * erfc));
    }
  }
  return check_result(al, gsl_sf_log_erfc(x));
}

static double amplgsl_sf_erf_Z(arglist *al) {
  double x = al->ra[0];
  double z = gsl_sf_erf_Z(x);
  if (al->derivs) {
    *al->derivs = -x * z;
    if (al->hes)
      *al->hes = -(z + x * *al->derivs);
  }
  return check_result(al, z);
}

static double amplgsl_sf_erf_Q(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = -gsl_sf_erf_Z(x);
    if (al->hes)
      *al->hes = -x * deriv;
  }
  return check_result(al, gsl_sf_erf_Q(x));
}

static double amplgsl_sf_hazard(arglist *al) {
  double x = al->ra[0];
  double hazard = gsl_sf_hazard(x);
  if (al->derivs) {
    *al->derivs = (hazard - x) * hazard;
    if (al->hes)
      *al->hes = hazard * (hazard * (2 * hazard - 3 * x) + (x * x - 1));
  }
  return check_result(al, hazard);
}

static double amplgsl_sf_expint_E1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -exp(-x) / x;
    if (al->hes)
      *al->hes = -*al->derivs * (1 / x + 1);
  }
  return check_result(al, gsl_sf_expint_E1(x));
}

static double amplgsl_sf_expint_E2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -gsl_sf_expint_E1(x);
    if (al->hes)
      *al->hes = exp(-x) / x;
  }
  return check_result(al, gsl_sf_expint_E2(x));
}

static double amplgsl_sf_expint_En(arglist *al) {
  int n = (int)al->ra[0];
  double x = al->ra[1];
  if (!check_int_arg(al, 0, "n"))
    return 0;
  if (al->derivs) {
    al->derivs[1] = n != 0 ?
        -gsl_sf_expint_En(n - 1, x) : -exp(-x) * (1 / x + 1) / x;
    if (al->hes) {
      if (n == 0)
        al->hes[2] = exp(-x) * (1 + 2 * (1 + 1 / x) / x) / x;
      else if (n == 1)
        al->hes[2] = exp(-x) * (1 / x + 1) / x;
      else
        al->hes[2] = gsl_sf_expint_En(n - 2, x);
    }
  }
  return check_result(al, gsl_sf_expint_En(n, x));
}

static double amplgsl_sf_expint_Ei(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(x) / x;
    if (al->hes)
      *al->hes = *al->derivs * (1 - 1 / x);
  }
  return check_result(al, gsl_sf_expint_Ei(x));
}

static double amplgsl_sf_Shi(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? sinh(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (cosh(x) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_Shi(x));
}

static double amplgsl_sf_Chi(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = cosh(x) / x;
    if (al->hes)
      *al->hes = (sinh(x) - *al->derivs) / x;
  }
  return check_result(al, gsl_sf_Chi(x));
}

static double amplgsl_sf_expint_3(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(-gsl_pow_3(x));
    if (al->hes)
      *al->hes = -3 * x * x * *al->derivs;
  }
  return check_result(al, gsl_sf_expint_3(x));
}

static double amplgsl_sf_Si(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? sin(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (cos(x) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_Si(x));
}

static double amplgsl_sf_Ci(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = cos(x) / x;
    if (al->hes)
      *al->hes = -(sin(x) + *al->derivs) / x;
  }
  return check_result(al, gsl_sf_Ci(x));
}

static double amplgsl_sf_atanint(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? atan(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (1 / (x * x + 1) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_atanint(x));
}

static double amplgsl_sf_fermi_dirac_m1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(x) / gsl_pow_2(exp(x) + 1);
    if (al->hes)
      *al->hes = -(exp(x) * (exp(x) - 1)) / gsl_pow_3(exp(x) + 1);
  }
  return check_result(al, gsl_sf_fermi_dirac_m1(x));
}

static double amplgsl_sf_fermi_dirac_0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_m1(x);
    if (al->hes)
      *al->hes = exp(x) / gsl_pow_2(exp(x) + 1);
  }
  return check_result(al, gsl_sf_fermi_dirac_0(x));
}

static double amplgsl_sf_fermi_dirac_1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_0(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_m1(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_1(x));
}

static double amplgsl_sf_fermi_dirac_2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_1(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_0(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_2(x));
}

static double amplgsl_sf_fermi_dirac_int(arglist *al) {
  int j = (int)al->ra[0];
  double x = al->ra[1];
  if (!check_int_arg(al, 0, "j"))
    return 0;
  if (al->derivs) {
    al->derivs[1] = gsl_sf_fermi_dirac_int(j - 1, x);
    if (al->hes)
      al->hes[2] = gsl_sf_fermi_dirac_int(j - 2, x);
  }
  return check_result(al, gsl_sf_fermi_dirac_int(j, x));
}

WRAP_CHECKED(gsl_sf_fermi_dirac_mhalf, ARGS1)
WRAP_CHECKED(gsl_sf_fermi_dirac_half, ARGS1)

static double amplgsl_sf_fermi_dirac_3half(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_half(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_mhalf(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_3half(x));
}

static double amplgsl_sf_fermi_dirac_inc_0(arglist *al) {
  double x = al->ra[0], b = al->ra[1];
  if (al->derivs) {
    double exp_x = exp(x), exp_b = exp(b);
    al->derivs[0] = exp_x / (exp_b + exp_x);
    al->derivs[1] = -al->derivs[0];
    if (al->hes) {
      al->hes[0] = al->hes[2] = al->derivs[0] * exp_b / (exp_b + exp_x);
      al->hes[1] = -al->hes[0];
    }
  }
  return check_result(al, gsl_sf_fermi_dirac_inc_0(x, b));
}

static double amplgsl_sf_gamma(arglist *al) {
  double x = al->ra[0];
  double gamma = gsl_sf_gamma(x);
  if (al->derivs) {
    double psi0 = gsl_sf_psi(x);
    *al->derivs = gamma * psi0;
    if (al->hes)
      *al->hes = *al->derivs * psi0 + gamma * gsl_sf_psi_1(x);
  }
  return check_result(al, gamma);
}

static double amplgsl_sf_lngamma(arglist *al) {
  double value = 0;
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x >= 0 || ceil(x) != x ? gsl_sf_psi(x) : GSL_NAN;
    if (al->hes)
      *al->hes = gsl_sf_psi_1(x);
  }
  CHECK_CALL(value, gsl_sf_lngamma_e(x, &result));
  return check_result(al, value);
}

static double amplgsl_sf_gammastar(arglist *al) {
  double x = al->ra[0];
  double gammastar = gsl_sf_gammastar(x);
  if (al->derivs) {
    double coef = (0.5 / x - log(x) + gsl_sf_psi(x));
    *al->derivs = coef * gammastar;
    if (al->hes) {
      *al->hes = coef * *al->derivs +
          (gsl_sf_psi_1(x) - (1 + 0.5 / x) / x) * gammastar;
    }
  }
  return check_result(al, gammastar);
}

static double amplgsl_sf_gammainv(arglist *al) {
  double x = al->ra[0];
  double gammainv = gsl_sf_gammainv(x);
  if (al->derivs) {
    if (x > 0 || ceil(x) != x) {
      double psi0 = gsl_sf_psi(x);
      *al->derivs = -gammainv * psi0;
      if (al->hes)
        *al->hes = -*al->derivs * psi0 - gammainv * gsl_sf_psi_1(x);
    } else {
      *al->derivs = pow(-1, -x) * gsl_sf_gamma(1 - x);
      if (al->hes)
        *al->hes = -2 * *al->derivs * gsl_sf_psi(1 - x);
    }
  }
  return check_result(al, gammainv);
}

WRAP_CHECKED(gsl_sf_poch, ARGS2)
WRAP_CHECKED(gsl_sf_lnpoch, ARGS2)
WRAP_CHECKED(gsl_sf_pochrel, ARGS2)

static double amplgsl_sf_gamma_inc(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs && check_const_arg(al, 0, "a")) {
    al->derivs[1] = x != 0 ? -exp(-x) * pow(x, a - 1) : GSL_NAN;
    if (al->hes)
      al->hes[2] = al->derivs[1] * (a - x - 1) / x;
  }
  return check_result(al, gsl_sf_gamma_inc(a, x));
}

WRAP_CHECKED(gsl_sf_gamma_inc_Q, ARGS2)
WRAP_CHECKED(gsl_sf_gamma_inc_P, ARGS2)

static double amplgsl_sf_beta(arglist *al) {
  double a = al->ra[0], b = al->ra[1];
  double beta = 0;
  CHECK_CALL(beta, gsl_sf_beta_e(a, b, &result));
  if (al->derivs) {
    double psi_a_plus_b = gsl_sf_psi(a + b);
    double da_coef = 0, db_coef = 0;
    int need_da = 1, need_db = 1;
    if (al->dig) {
      need_da = !al->dig[0];
      need_db = !al->dig[1];
    }
    if (need_da) {
      da_coef = gsl_sf_psi(a) - psi_a_plus_b;
      al->derivs[0] = beta * da_coef;
    }
    if (need_db) {
      db_coef = gsl_sf_psi(b) - psi_a_plus_b;
      al->derivs[1] = beta * db_coef;
    }
    if (al->hes) {
      double psi1_a_plus_b = gsl_sf_psi_1(a + b);
      if (need_da) {
        al->hes[0] = al->derivs[0] * da_coef +
            beta * (gsl_sf_psi_1(a) - psi1_a_plus_b);
        if (need_db)
          al->hes[1] = al->derivs[0] * db_coef - beta * psi1_a_plus_b;
      }
      if (need_db) {
       al->hes[2] = al->derivs[1] * db_coef +
           beta * (gsl_sf_psi_1(b) - psi1_a_plus_b);
      }
    }
  }
  return check_result(al, beta);
}

static double amplgsl_sf_lnbeta(arglist *al) {
  double a = al->ra[0], b = al->ra[1];
  if (al->derivs) {
    double psi_a_plus_b = gsl_sf_psi(a + b);
    int need_da = 1, need_db = 1;
    if (al->dig) {
      need_da = !al->dig[0];
      need_db = !al->dig[1];
    }
    if (need_da)
      al->derivs[0] = gsl_sf_psi(a) - psi_a_plus_b;
    if (need_db)
      al->derivs[1] = gsl_sf_psi(b) - psi_a_plus_b;
    if (al->hes) {
      double psi1_a_plus_b = gsl_sf_psi_1(a + b);
      if (need_da) {
        al->hes[0] = gsl_sf_psi_1(a) - psi1_a_plus_b;
        if (need_db)
          al->hes[1] = -psi1_a_plus_b;
      }
      if (need_db)
       al->hes[2] = gsl_sf_psi_1(b) - psi1_a_plus_b;
    }
  }
  return check_result(al, gsl_sf_lnbeta(a, b));
}

WRAP_CHECKED(gsl_sf_beta_inc, ARGS3)

static double amplgsl_sf_gegenpoly_1(arglist *al) {
  double lambda = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    al->derivs[0] = 2 * x;
    /* For unclear reason gsl_sf_gegenpoly_1(0, x) returns 2 * x. */
    al->derivs[1] = lambda != 0 ? 2 * lambda : 2;
    if (al->hes) {
      al->hes[0] = al->hes[2] = 0;
      al->hes[1] = 2;
    }
  }
  return check_result(al, gsl_sf_gegenpoly_1(lambda, x));
}

static double amplgsl_sf_gegenpoly_2(arglist *al) {
  double lambda = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    double coef1 = (1 + 2 * lambda) * x;
    double coef2 = 4 * lambda * (lambda + 1);
    al->derivs[0] = 2 * coef1 * x - 1;
    /* For unclear reason gsl_sf_gegenpoly_2(0, x) returns 2 * x^2 - 1. */
    al->derivs[1] = lambda != 0 ? coef2 * x : 4 * x;
    if (al->hes) {
      al->hes[0] = 4 * x * x;
      al->hes[1] = 4 * coef1;
      al->hes[2] = lambda != 0 ? coef2 : 4;
    }
  }
  return check_result(al, gsl_sf_gegenpoly_2(lambda, x));
}

static double amplgsl_sf_gegenpoly_3(arglist *al) {
  double lambda = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    double x2 = x * x;
    al->derivs[0] =
        x * (4 * (2.0 / 3.0 + lambda * (lambda + 2)) * x2 -
        2 * (2 * lambda + 1));
    /* For unclear reason gsl_sf_gegenpoly_3(0, x) returns
       x * (-2.0 + 4.0 / 3.0 * x * x). */
    al->derivs[1] = lambda != 0 ?
        2 * lambda * (lambda + 1) * (2 * (2 + lambda) * x2 - 1) :
        4 * x2 - 2;
    if (al->hes) {
      al->hes[0] = 4 * x * (2 * (lambda + 1) * x2 - 1);
      al->hes[1] = 2 * (x2 * (6 * lambda * lambda + 4) +
          2 * lambda * (6 * x2 - 1) - 1);
      al->hes[2] = lambda != 0 ?
          8 * lambda * (lambda + 1) * (lambda + 2) * x : 8 * x;
    }
  }
  return check_result(al, gsl_sf_gegenpoly_3(lambda, x));
}

static double amplgsl_sf_gegenpoly_n(arglist *al) {
  int n = (int)al->ra[0];
  double lambda = al->ra[1], x = al->ra[2];
  if (!check_int_arg(al, 0, "n"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al, gsl_sf_gegenpoly_n(n, lambda, x));
}

static double amplgsl_sf_hyperg_0F1(arglist *al) {
  double c = 0, x = 0;
  if (!check_args(al))
    return 0;
  c = al->ra[0];
  x = al->ra[1];
  if (al->derivs && check_const_arg(al, 0, "c")) {
    al->derivs[1] = gsl_sf_hyperg_0F1(c + 1, x) / c;
    if (al->hes)
      al->hes[2] = gsl_sf_hyperg_0F1(c + 2, x) / (c * (c + 1));
  }
  return check_result(al, gsl_sf_hyperg_0F1(c, x));
}

static double amplgsl_sf_hyperg_1F1_int(arglist *al) {
  int m = 0, n = 0;
  double x = 0;
  if (!check_args(al) || !check_int_arg(al, 0, "m") ||
      !check_int_arg(al, 1, "n")) {
    return 0;
  }
  m = (int)al->ra[0];
  n = (int)al->ra[1];
  x = al->ra[2];
  if (al->derivs) {
    /* If n is an integer <= 0, then 1F1(m; n; x) is undefined.
       See http://mathworld.wolfram.com/
       ConfluentHypergeometricFunctionoftheFirstKind.html */
    al->derivs[2] = n > 0 ?
        m * gsl_sf_hyperg_1F1_int(m + 1, n + 1, x) / n : GSL_NAN;
    if (al->hes) {
      al->hes[5] =
          m * (m + 1) * gsl_sf_hyperg_1F1_int(m + 2, n + 2, x) / (n * (n + 1));
    }
  }
  return check_result(al, gsl_sf_hyperg_1F1_int(m, n, x));
}

WRAP_CHECKED(gsl_sf_hyperg_1F1, ARGS3)

static double amplgsl_sf_hyperg_U_int(arglist *al) {
  double value = 0;
  if (!check_args(al) || !check_int_arg(al, 0, "m") ||
      !check_int_arg(al, 1, "n")) {
    return 0;
  }
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  CHECK_CALL(value, gsl_sf_hyperg_U_int_e(
    (int)al->ra[0], (int)al->ra[1], al->ra[2], &result));
  return check_result(al, value);
}

WRAP_CHECKED(gsl_sf_hyperg_U, ARGS3)
WRAP_CHECKED(gsl_sf_hyperg_2F1, ARGS4)
WRAP_CHECKED(gsl_sf_hyperg_2F1_conj, ARGS4)
WRAP_CHECKED(gsl_sf_hyperg_2F1_renorm, ARGS4)
WRAP_CHECKED(gsl_sf_hyperg_2F1_conj_renorm, ARGS4)
WRAP_CHECKED(gsl_sf_hyperg_2F0, ARGS3)

static double amplgsl_sf_laguerre_1(arglist *al) {
  double a = 0, x = 0;
  if (!check_args(al))
    return 0;
  a = al->ra[0];
  x = al->ra[1];
  if (al->derivs) {
    al->derivs[0] =  1;
    al->derivs[1] = -1;
    if (al->hes)
      al->hes[0] = al->hes[1] = al->hes[2] = 0;
  }
  return check_result(al, gsl_sf_laguerre_1(a, x));
}

static double amplgsl_sf_laguerre_2(arglist *al) {
  double a = 0, x = 0;
  if (!check_args(al))
    return 0;
  a = al->ra[0];
  x = al->ra[1];
  if (al->derivs) {
    al->derivs[0] =  a - x + 1.5;
    al->derivs[1] = -a + x - 2;
    if (al->hes) {
      al->hes[0] = al->hes[2] = 1;
      al->hes[1] = -1;
    }
  }
  return check_result(al, gsl_sf_laguerre_2(a, x));
}

static double amplgsl_sf_laguerre_3(arglist *al) {
  double a = 0, x = 0;
  if (!check_args(al))
    return 0;
  a = al->ra[0];
  x = al->ra[1];
  if (al->derivs) {
    al->derivs[0] = (11 + 3 * a * (a - 2 * (x - 2)) + 3 * x * (x - 5)) / 6;
    al->derivs[1] = (-0.5 * (a + 2) + x) * (a + 3) - 0.5 * x * x;
    if (al->hes) {
      al->hes[0] =  a - x + 2;
      al->hes[1] = -a + x - 2.5;
      al->hes[2] =  a - x + 3;
    }
  }
  return check_result(al, gsl_sf_laguerre_3(a, x));
}

static double amplgsl_sf_laguerre_n(arglist *al) {
  if (!check_args(al) || !check_int_arg(al, 0, "n"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_laguerre_n((int)al->ra[0], al->ra[1], al->ra[2]));
}

static double amplgsl_sf_lambert_W0(arglist *al) {
  double x = al->ra[0];
  double value = 0;
  CHECK_CALL(value, gsl_sf_lambert_W0_e(x, &result));
  if (al->derivs) {
    if (x < -1 / M_E)
      *al->derivs = GSL_NAN;
    else
      *al->derivs = x != 0 ? value / (x * (value + 1)) : 1;
    if (al->hes)
      *al->hes = -*al->derivs * *al->derivs * (value + 2) / (value + 1);
  }
  return check_result(al, value);
}

static double amplgsl_sf_lambert_Wm1(arglist *al) {
  double x = al->ra[0];
  double value = 0;
  CHECK_CALL(value, gsl_sf_lambert_Wm1_e(x, &result));
  if (al->derivs) {
    if (x < -1 / M_E)
      *al->derivs = GSL_NAN;
    else
      *al->derivs = value / (x * (value + 1));
    if (al->hes)
      *al->hes = -*al->derivs * *al->derivs * (value + 2) / (value + 1);
  }
  return check_result(al, value);
}

static double amplgsl_sf_legendre_P1(arglist *al) {
  if (al->derivs) {
    *al->derivs = 1;
    if (al->hes)
      *al->hes = 0;
  }
  return check_result(al, gsl_sf_legendre_P1(al->ra[0]));
}

static double amplgsl_sf_legendre_P2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 3 * x;
    if (al->hes)
      *al->hes = 3;
  }
  return check_result(al, gsl_sf_legendre_P2(x));
}

static double amplgsl_sf_legendre_P3(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 7.5 * x * x - 1.5;
    if (al->hes)
      *al->hes = 15 * x;
  }
  return check_result(al, gsl_sf_legendre_P3(x));
}

static double amplgsl_sf_legendre_Pl(arglist *al) {
  int el = 0;
  double x = 0, pl = 0;
  if (!check_int_arg(al, 0, "l"))
    return 0;
  el = (int)al->ra[0];
  x = al->ra[1];
  pl = gsl_sf_legendre_Pl(el, x);
  if (al->derivs) {
    if (fabs(x) != 1) {
      double pl_plus_1 = gsl_sf_legendre_Pl(el + 1, x);
      double coef = (el + 1) / (x * x - 1);
      al->derivs[1] = -coef * (x * pl - pl_plus_1);
      if (al->hes) {
        al->hes[2] =
            coef * ((x * x * (el + 2) + 1) * pl - (2 * el + 5) * x * pl_plus_1 +
            (el + 2) * gsl_sf_legendre_Pl(el + 2, x)) / (x * x - 1);
      }
    } else {
      double coef = 0.5 * el * (el + 1);
      al->derivs[1] = pow(x, el + 1) * coef;
      if (al->hes)
        al->hes[2] = pow(x, el) * 0.25 * coef * (el * el + el - 2);
    }
  }
  return check_result(al, pl);
}

static double amplgsl_sf_legendre_Q0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 1 / (1 - x * x);
    if (al->hes)
      *al->hes = 2 * x * *al->derivs * *al->derivs;
  }
  return check_result(al, gsl_sf_legendre_Q0(x));
}

static double amplgsl_sf_legendre_Q1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double coef = 1 / (1 - x * x);
    *al->derivs = coef * x + 0.5 * (log(1 + x) - log(fabs(1 - x)));
    if (al->hes)
      *al->hes = 2 * coef * coef;
  }
  return check_result(al, gsl_sf_legendre_Q1(x));
}

static double amplgsl_sf_legendre_Ql(arglist *al) {
  int el = 0;
  double x = 0, ql = 0;
  if (!check_int_arg(al, 0, "l"))
    return 0;
  el = (int)al->ra[0];
  x = al->ra[1];
  ql = gsl_sf_legendre_Ql(el, x);
  if (al->derivs) {
    double coef = (el + 1) / (x * x - 1);
    double ql_plus_1 = gsl_sf_legendre_Ql(el + 1, x);
    al->derivs[1] = coef * (ql_plus_1 - x * ql);
    if (al->hes) {
      al->hes[2] =
          coef * ((x * x * (el + 2) + 1) * ql - (2 * el + 5) * x * ql_plus_1 +
          (el + 2) * gsl_sf_legendre_Ql(el + 2, x)) / (x * x - 1);
    }
  }
  return check_result(al, ql);
}

static double amplgsl_sf_legendre_Plm(arglist *al) {
  if (!check_int_arg(al, 0, "l") || !check_int_arg(al, 1, "m"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_legendre_Plm((int)al->ra[0], (int)al->ra[1], al->ra[2]));
}

static double amplgsl_sf_legendre_sphPlm(arglist *al) {
  if (!check_int_arg(al, 0, "l") || !check_int_arg(al, 1, "m"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_legendre_sphPlm((int)al->ra[0], (int)al->ra[1], al->ra[2]));
}

WRAP_CHECKED(gsl_sf_conicalP_half, ARGS2)
WRAP_CHECKED(gsl_sf_conicalP_mhalf, ARGS2)
WRAP_CHECKED(gsl_sf_conicalP_0, ARGS2)
WRAP_CHECKED(gsl_sf_conicalP_1, ARGS2)

static double amplgsl_sf_conicalP_sph_reg(arglist *al) {
  double value = 0;
  if (!check_int_arg(al, 0, "m"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  CHECK_CALL(value, gsl_sf_conicalP_sph_reg_e(
    (int)al->ra[0], al->ra[1], al->ra[2], &result));
  return check_result(al, value);
}

static double amplgsl_sf_conicalP_cyl_reg(arglist *al) {
  double value = 0;
  if (!check_int_arg(al, 0, "m"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  CHECK_CALL(value, gsl_sf_conicalP_cyl_reg_e(
      (int)al->ra[0], al->ra[1], al->ra[2], &result));
  return check_result(al, value);
}

WRAP_CHECKED(gsl_sf_legendre_H3d_0, ARGS2)
WRAP_CHECKED(gsl_sf_legendre_H3d_1, ARGS2)

static double amplgsl_sf_legendre_H3d(arglist *al) {
  double value = 0;
  if (!check_int_arg(al, 0, "l"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  CHECK_CALL(value, gsl_sf_legendre_H3d_e(
      (int)al->ra[0], al->ra[1], al->ra[2], &result));
  return check_result(al, value);
}

static double amplgsl_sf_log(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = 1 / x;
    if (al->hes)
      *al->hes = -deriv * deriv;
  }
  return check_result(al, gsl_sf_log(x));
}

static double amplgsl_sf_log_abs(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = 1 / x;
    if (al->hes)
      *al->hes = -deriv * deriv;
  }
  return check_result(al, gsl_sf_log_abs(x));
}

static double amplgsl_sf_log_1plusx(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = 1 / (1 + x);
    if (al->hes)
      *al->hes = -deriv * deriv;
  }
  return check_result(al, gsl_sf_log_1plusx(x));
}

static double amplgsl_sf_log_1plusx_mx(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double sub = 1 / (1 + x);
    *al->derivs = sub - 1;
    if (al->hes)
      *al->hes = -sub * sub;
  }
  return check_result(al, gsl_sf_log_1plusx_mx(x));
}

static double amplgsl_sf_mathieu_a(arglist *al) {
  int n = 0;
  double q = 0;
  gsl_sf_result result = {0, 0};
  if (!check_int_arg(al, 0, "n"))
    return 0;
  n = (int)al->ra[0];
  q = al->ra[1];
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_mathieu_a_e(n, q, &result) ? GSL_NAN : result.val);
}

static double amplgsl_sf_mathieu_b(arglist *al) {
  int n = 0;
  double q = 0;
  gsl_sf_result result = {0, 0};
  if (!check_int_arg(al, 0, "n"))
    return 0;
  n = (int)al->ra[0];
  q = al->ra[1];
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_mathieu_b_e(n, q, &result) ? GSL_NAN : result.val);
}

static double amplgsl_sf_mathieu_ce(arglist *al) {
  int n = 0;
  double q = 0, x = 0;
  gsl_sf_result result = {0, 0};
  if (!check_int_arg(al, 0, "n") || !check_args(al))
    return 0;
  n = (int)al->ra[0];
  q = al->ra[1];
  x = al->ra[2];
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_mathieu_ce_e(n, q, x, &result) ? GSL_NAN : result.val);
}

static double amplgsl_sf_mathieu_se(arglist *al) {
  int n = 0;
  double q = 0, x = 0;
  gsl_sf_result result = {0, 0};
  if (!check_int_arg(al, 0, "n"))
    return 0;
  n = (int)al->ra[0];
  q = al->ra[1];
  x = al->ra[2];
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_mathieu_se_e(n, q, x, &result) ? GSL_NAN : result.val);
}

static double amplgsl_sf_mathieu_Mc(arglist *al) {
  int j = 0, n = 0;
  double q = 0, x = 0;
  gsl_sf_result result = {0, 0};
  if (!check_int_arg(al, 0, "j") || !check_int_arg(al, 1, "n"))
    return 0;
  j = (int)al->ra[0];
  n = (int)al->ra[1];
  q = al->ra[2];
  x = al->ra[3];
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_mathieu_Mc_e(j, n, q, x, &result) ? GSL_NAN : result.val);
}

static double amplgsl_sf_mathieu_Ms(arglist *al) {
  int j = 0, n = 0;
  double q = 0, x = 0;
  gsl_sf_result result = {0, 0};
  if (!check_int_arg(al, 0, "j") || !check_int_arg(al, 1, "n"))
    return 0;
  j = (int)al->ra[0];
  n = (int)al->ra[1];
  q = al->ra[2];
  x = al->ra[3];
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_sf_mathieu_Ms_e(j, n, q, x, &result) ? GSL_NAN : result.val);
}

static double amplgsl_sf_pow_int(arglist *al) {
  double value = 0;
  double x = 0;
  int n = 0;
  if (!check_int_arg(al, 1, "n"))
    return 0;
  x = al->ra[0];
  n = (int)al->ra[1];
  if (al->derivs) {
    *al->derivs = n != 0 ? n * gsl_sf_pow_int(x, n - 1) : 0;
    if (al->hes)
      *al->hes = n != 0 && n != 1 ? n * (n - 1) * gsl_sf_pow_int(x, n - 2) : 0;
  }
  CHECK_CALL(value, gsl_sf_pow_int_e(x, n, &result));
  return check_result(al, value);
}

static double amplgsl_sf_psi_int(arglist *al) {
  if (!check_int_arg(al, 0, "n"))
    return 0;
  return check_result(al, gsl_sf_psi_int((int)al->ra[0]));
}

static double amplgsl_sf_psi(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x >= 0 || ceil(x) != x ? gsl_sf_psi_1(x) : GSL_NAN;
    if (al->hes)
      *al->hes = gsl_sf_psi_n(2, x);
  }
  return check_result(al, gsl_sf_psi(x));
}

WRAP_CHECKED(gsl_sf_psi_1piy, ARGS1)

static double amplgsl_sf_psi_1_int(arglist *al) {
  if (!check_int_arg(al, 0, "n"))
    return 0;
  return check_result(al, gsl_sf_psi_1_int((int)al->ra[0]));
}

static double amplgsl_sf_psi_1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_psi_n(2, x);
    if (al->hes)
      *al->hes = gsl_sf_psi_n(3, x);
  }
  return check_result(al, gsl_sf_psi_1(x));
}

static double amplgsl_sf_psi_n(arglist *al) {
  int n = 0;
  double x = 0;
  if (!check_int_arg(al, 0, "n"))
    return 0;
  n = (int)al->ra[0];
  x = al->ra[1];
  if (al->derivs) {
    al->derivs[1] = x >= 0 || ceil(x) != x ? gsl_sf_psi_n(n + 1, x) : GSL_NAN;
    if (al->hes)
      al->hes[2] = gsl_sf_psi_n(n + 2, x);
  }
  return check_result(al, gsl_sf_psi_n(n, x));
}

WRAP_CHECKED(gsl_sf_synchrotron_1, ARGS1)
WRAP_CHECKED(gsl_sf_synchrotron_2, ARGS1)

static double amplgsl_sf_transport_2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    if (x != 0) {
      double exp_x = exp(x);
      double coef = exp_x * x / gsl_pow_2(exp_x - 1);
      *al->derivs = coef * x;
      if (al->hes)
        *al->hes = -coef * (exp_x * (x - 2) + x + 2) / (exp_x - 1);
    } else {
      *al->derivs = 1;
      if (al->hes)
        *al->hes = 0;
    }
  }
  return check_result(al, gsl_sf_transport_2(x));
}

static double amplgsl_sf_transport_3(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    if (x != 0) {
      double exp_x = exp(x);
      double coef = exp_x * x * x / gsl_pow_2(exp_x - 1);
      *al->derivs = coef * x;
      if (al->hes)
        *al->hes = -coef * (exp_x * (x - 3) + x + 3) / (exp_x - 1);
    } else {
      *al->derivs = 0;
      if (al->hes)
        *al->hes = 1;
    }
  }
  return check_result(al, gsl_sf_transport_3(x));
}

static double amplgsl_sf_transport_4(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    if (x != 0) {
      double exp_x = exp(x);
      double coef = exp_x * x * x * x / gsl_pow_2(exp_x - 1);
      *al->derivs = coef * x;
      if (al->hes)
        *al->hes = -coef * (exp_x * (x - 4) + x + 4) / (exp_x - 1);
    } else {
      *al->derivs = 0;
      if (al->hes)
        *al->hes = 0;
    }
  }
  return check_result(al, gsl_sf_transport_4(x));
}

static double amplgsl_sf_transport_5(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    if (x != 0) {
      double exp_x = exp(x);
      double x2 = x * x;
      double coef = exp_x * x2 * x2 / gsl_pow_2(exp_x - 1);
      *al->derivs = coef * x;
      if (al->hes)
        *al->hes = -coef * (exp_x * (x - 5) + x + 5) / (exp_x - 1);
    } else {
      *al->derivs = 0;
      if (al->hes)
        *al->hes = 0;
    }
  }
  return check_result(al, gsl_sf_transport_5(x));
}


static double amplgsl_sf_sinc(arglist* al) {
  double x = al->ra[0];
  double sinc = gsl_sf_sinc(x);
  if (x != 0)
  {
    double pix = M_PI * x;
    double cospi = cos(pix);
    double xsquare = gsl_sf_pow_int(x, 2);
    if (al->derivs) {
      *al->derivs = ((cospi * M_PI) - sinc) / x;
      if (al->hes) {
        *al->hes = -(sinc * gsl_sf_pow_int(M_PI,2)) + 2 * (sinc / xsquare) 
          - 2 * (cospi) / xsquare;
      }
    }
  }
  else
  {
    if (al->derivs) {
      *al->derivs = 0;
      if (al->hes) {
        *al->hes = -1 / 3;
      }
    }
  }
  return check_result(al, sinc);
}

static double amplgsl_sf_zeta_int(arglist *al) {
  
  if (!check_int_arg(al, 0, "n"))
    return 0;
  return check_result(al, gsl_sf_zeta_int((int)al->ra[0]));
}

WRAP_CHECKED(gsl_sf_zeta, ARGS1)

static double amplgsl_sf_zetam1_int(arglist *al) {
  if (!check_int_arg(al, 0, "n"))
    return 0;
  return check_result(al, gsl_sf_zetam1_int((int)al->ra[0]));
}

WRAP_CHECKED(gsl_sf_zetam1, ARGS1)
WRAP_CHECKED(gsl_sf_hzeta, ARGS2)

static double amplgsl_sf_eta_int(arglist *al) {
  if (!check_int_arg(al, 0, "n"))
    return 0;
  return check_result(al, gsl_sf_eta_int((int)al->ra[0]));
}

WRAP_CHECKED(gsl_sf_eta, ARGS1)

static gsl_rng *rng;

static void free_rng(void *data) {
  UNUSED(data);
  if (rng) {
	gsl_rng_free(rng);
	rng = 0;
	}
}

#ifdef addrandinit
 static void
rng_init(void *v, unsigned long x)
{
	UNUSED(v);
	gsl_rng_default_seed = x;
	if (rng)
		gsl_rng_free(rng);
	rng = gsl_rng_alloc(gsl_rng_env_setup());
	}
#endif

WRAP(gsl_ran_gaussian, RNG_ARGS1)

static double amplgsl_ran_gaussian_pdf(arglist *al) {
  double x = al->ra[0], sigma = al->ra[1];
  double pdf = gsl_ran_gaussian_pdf(x, sigma);
  if (al->derivs) {
    double sigma2 = sigma * sigma, sigma3 = sigma2 * sigma, x2 = x * x;
    al->derivs[0] = -x * pdf / sigma2;
    al->derivs[1] = (x2 - sigma2) * pdf / sigma3;
    if (al->hes) {
      al->hes[0] = al->derivs[1] / sigma;
      al->hes[1] = (x2 - 3 * sigma2) * al->derivs[0] / sigma3;
      al->hes[2] = (x2 * (x2 / sigma2 - 5) / sigma2 + 2) * pdf / sigma2;
    }
  }
  return check_result(al, pdf);
}

WRAP(gsl_ran_gaussian_ziggurat, RNG_ARGS1)
WRAP(gsl_ran_gaussian_ratio_method, RNG_ARGS1)
WRAP(gsl_ran_ugaussian, rng)

static double amplgsl_ran_ugaussian_pdf(arglist *al) {
  double x = al->ra[0];
  double pdf = gsl_ran_ugaussian_pdf(x);
  if (al->derivs) {
    *al->derivs = -x * pdf;
    if (al->hes)
      *al->hes = (x * x - 1) * pdf;
  }
  return check_result(al, pdf);
}

WRAP(gsl_ran_ugaussian_ratio_method, rng)

static double amplgsl_cdf_gaussian_P(arglist *al) {
  double x = al->ra[0], sigma = al->ra[1];
  if (al->derivs) {
    double pdf = gsl_ran_gaussian_pdf(x, sigma);
    double sign = sigma >= 0 ? 1 : -1;
    al->derivs[0] = sign * pdf;
    al->derivs[1] = sign * -x * pdf / sigma;
    if (al->hes) {
      double x2 = x * x;
      double sigma2 = sigma * sigma, sigma3 = sigma2 * sigma;
      al->hes[0] = al->derivs[1] / sigma;
      al->hes[1] = (x2 - sigma2) * al->derivs[0] / sigma3;
      al->hes[2] = (x2 - 2 * sigma2) * al->derivs[1] / sigma3;
    }
  }
  return check_result(al, gsl_cdf_gaussian_P(x, sigma));
}

WRAP(gsl_cdf_gaussian_Q, ARGS2)
WRAP(gsl_cdf_gaussian_Pinv, ARGS2)
WRAP(gsl_cdf_gaussian_Qinv, ARGS2)

static double amplgsl_cdf_ugaussian_P(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double pdf = gsl_ran_ugaussian_pdf(x);
    *al->derivs = pdf;
    if (al->hes)
      *al->hes = -x * pdf;
  }
  return check_result(al, gsl_cdf_ugaussian_P(x));
}

WRAP(gsl_cdf_ugaussian_Q, ARGS1)

static double amplgsl_cdf_ugaussian_Pinv(arglist *al) {
  double x = al->ra[0];
  double pinv = gsl_cdf_ugaussian_Pinv(x);
  if (al->derivs) {
    double part = 2 * M_PI * exp(pinv * pinv);
    *al->derivs = sqrt(part);
    if (al->hes)
      *al->hes = part * pinv;
  }
  return check_result(al, pinv);
}

WRAP(gsl_cdf_ugaussian_Qinv, ARGS1)
WRAP(gsl_ran_gaussian_tail, RNG_ARGS2)
WRAP(gsl_ran_gaussian_tail_pdf, ARGS3)
WRAP(gsl_ran_ugaussian_tail, RNG_ARGS1)
WRAP(gsl_ran_ugaussian_tail_pdf, ARGS2)

WRAP(gsl_ran_exponential, RNG_ARGS1)

static double amplgsl_ran_exponential_pdf(arglist *al) {
  double x = al->ra[0], mu = al->ra[1];
  double pdf = gsl_ran_exponential_pdf(x, mu);
  if (al->derivs) {
    double mu2 = mu * mu;
    al->derivs[0] = x < 0 ? 0 : -pdf / mu;
    al->derivs[1] = x < 0 ? 0 : (x - mu) * pdf / mu2;
    if (al->hes) {
      al->hes[0] = x < 0 ? 0 : -al->derivs[0] / mu;
      al->hes[1] = x < 0 ? 0 : (2 * mu - x) * al->hes[0] / mu;
      al->hes[2] = x < 0 ? 0 : ((x * x / mu - 4 * x) / mu + 2) * al->hes[0];
    }
  }
  return check_result(al, pdf);
}

WRAP(gsl_cdf_exponential_P, ARGS2)
WRAP(gsl_cdf_exponential_Q, ARGS2)
WRAP(gsl_cdf_exponential_Pinv, ARGS2)
WRAP(gsl_cdf_exponential_Qinv, ARGS2)

WRAP(gsl_ran_laplace, RNG_ARGS1)

static double amplgsl_ran_laplace_pdf(arglist *al) {
  double x = al->ra[0], a = al->ra[1];
  double pdf = gsl_ran_laplace_pdf(x, a);
  if (al->derivs) {
    double dx = pdf / a;
    if (!al->dig || !al->dig[0]) {
      al->derivs[0] = -x * dx / fabs(x);
      if (al->hes) {
        al->hes[0] = dx / a;
        al->hes[1] = mul_by_sign(x, 2 - fabs(x) / a) * al->hes[0];
      }
    }
    al->derivs[1] = (fabs(x) - a) * dx / a;
    if (al->hes)
      al->hes[2] = ((x * x / a - 4 * fabs(x)) / a + 2) * dx / a;
  }
  return check_result(al, pdf);
}

WRAP(gsl_cdf_laplace_P, ARGS2)
WRAP(gsl_cdf_laplace_Q, ARGS2)
WRAP(gsl_cdf_laplace_Pinv, ARGS2)
WRAP(gsl_cdf_laplace_Qinv, ARGS2)

WRAP(gsl_ran_exppow, RNG_ARGS2)
WRAP(gsl_ran_exppow_pdf, ARGS3)
WRAP(gsl_cdf_exppow_P, ARGS3)
WRAP(gsl_cdf_exppow_Q, ARGS3)

WRAP(gsl_ran_cauchy, RNG_ARGS1)
WRAP(gsl_ran_cauchy_pdf, ARGS2)
WRAP(gsl_cdf_cauchy_P, ARGS2)
WRAP(gsl_cdf_cauchy_Q, ARGS2)
WRAP(gsl_cdf_cauchy_Pinv, ARGS2)
WRAP(gsl_cdf_cauchy_Qinv, ARGS2)

WRAP(gsl_ran_rayleigh, RNG_ARGS1)
WRAP(gsl_ran_rayleigh_pdf, ARGS2)
WRAP(gsl_cdf_rayleigh_P, ARGS2)
WRAP(gsl_cdf_rayleigh_Q, ARGS2)
WRAP(gsl_cdf_rayleigh_Pinv, ARGS2)
WRAP(gsl_cdf_rayleigh_Qinv, ARGS2)

WRAP(gsl_ran_rayleigh_tail, RNG_ARGS2)
WRAP(gsl_ran_rayleigh_tail_pdf, ARGS3)

WRAP(gsl_ran_landau, rng)
WRAP(gsl_ran_landau_pdf, ARGS1)

WRAP(gsl_ran_levy, RNG_ARGS2)
WRAP(gsl_ran_levy_skew, RNG_ARGS3)

WRAP(gsl_ran_gamma, RNG_ARGS2)
WRAP(gsl_ran_gamma_knuth, RNG_ARGS2)
WRAP(gsl_ran_gamma_pdf, ARGS3)

static double amplgsl_cdf_gamma_P(arglist *al) {
  double x = al->ra[0], a = al->ra[1], b = al->ra[2];
  if (al->derivs && check_const_arg(al, 1, "a") &&
      check_const_arg(al, 2, "b") && !al->dig[0]) {
    *al->derivs = pow(x / b, a) / (exp(x / b) * x * gsl_sf_gamma(a));
    if (al->hes)
      *al->hes = *al->derivs * ((a - 1) / x - 1 / b);
  }
  return check_result(al, gsl_cdf_gamma_P(x, a, b));
}

WRAP(gsl_cdf_gamma_Q, ARGS3)
WRAP(gsl_cdf_gamma_Pinv, ARGS3)
WRAP(gsl_cdf_gamma_Qinv, ARGS3)

WRAP(gsl_ran_flat, RNG_ARGS2)
WRAP(gsl_ran_flat_pdf, ARGS3)
WRAP(gsl_cdf_flat_P, ARGS3)
WRAP(gsl_cdf_flat_Q, ARGS3)
WRAP(gsl_cdf_flat_Pinv, ARGS3)
WRAP(gsl_cdf_flat_Qinv, ARGS3)

WRAP(gsl_ran_lognormal, RNG_ARGS2)
WRAP(gsl_ran_lognormal_pdf, ARGS3)
WRAP(gsl_cdf_lognormal_P, ARGS3)
WRAP(gsl_cdf_lognormal_Q, ARGS3)
WRAP(gsl_cdf_lognormal_Pinv, ARGS3)
WRAP(gsl_cdf_lognormal_Qinv, ARGS3)

WRAP(gsl_ran_chisq, RNG_ARGS1)
WRAP(gsl_ran_chisq_pdf, ARGS2)
WRAP(gsl_cdf_chisq_P, ARGS2)
WRAP(gsl_cdf_chisq_Q, ARGS2)
WRAP(gsl_cdf_chisq_Pinv, ARGS2)
WRAP(gsl_cdf_chisq_Qinv, ARGS2)

WRAP(gsl_ran_fdist, RNG_ARGS2)
WRAP(gsl_ran_fdist_pdf, ARGS3)
WRAP(gsl_cdf_fdist_P, ARGS3)
WRAP(gsl_cdf_fdist_Q, ARGS3)
WRAP(gsl_cdf_fdist_Pinv, ARGS3)
WRAP(gsl_cdf_fdist_Qinv, ARGS3)

WRAP(gsl_ran_tdist, RNG_ARGS1)
WRAP(gsl_ran_tdist_pdf, ARGS2)
WRAP(gsl_cdf_tdist_P, ARGS2)
WRAP(gsl_cdf_tdist_Q, ARGS2)
WRAP(gsl_cdf_tdist_Pinv, ARGS2)
WRAP(gsl_cdf_tdist_Qinv, ARGS2)

WRAP(gsl_ran_beta, RNG_ARGS2)
WRAP(gsl_ran_beta_pdf, ARGS3)
WRAP(gsl_cdf_beta_P, ARGS3)
WRAP(gsl_cdf_beta_Q, ARGS3)
WRAP(gsl_cdf_beta_Pinv, ARGS3)
WRAP(gsl_cdf_beta_Qinv, ARGS3)

WRAP(gsl_ran_logistic, RNG_ARGS1)
WRAP(gsl_ran_logistic_pdf, ARGS2)
WRAP(gsl_cdf_logistic_P, ARGS2)
WRAP(gsl_cdf_logistic_Q, ARGS2)
WRAP(gsl_cdf_logistic_Pinv, ARGS2)
WRAP(gsl_cdf_logistic_Qinv, ARGS2)

WRAP(gsl_ran_pareto, RNG_ARGS2)
WRAP(gsl_ran_pareto_pdf, ARGS3)
WRAP(gsl_cdf_pareto_P, ARGS3)
WRAP(gsl_cdf_pareto_Q, ARGS3)
WRAP(gsl_cdf_pareto_Pinv, ARGS3)
WRAP(gsl_cdf_pareto_Qinv, ARGS3)

WRAP(gsl_ran_weibull, RNG_ARGS2)
WRAP(gsl_ran_weibull_pdf, ARGS3)
WRAP(gsl_cdf_weibull_P, ARGS3)
WRAP(gsl_cdf_weibull_Q, ARGS3)
WRAP(gsl_cdf_weibull_Pinv, ARGS3)
WRAP(gsl_cdf_weibull_Qinv, ARGS3)

WRAP(gsl_ran_gumbel1, RNG_ARGS2)
WRAP(gsl_ran_gumbel1_pdf, ARGS3)
WRAP(gsl_cdf_gumbel1_P, ARGS3)
WRAP(gsl_cdf_gumbel1_Q, ARGS3)
WRAP(gsl_cdf_gumbel1_Pinv, ARGS3)
WRAP(gsl_cdf_gumbel1_Qinv, ARGS3)

WRAP(gsl_ran_gumbel2, RNG_ARGS2)
WRAP(gsl_ran_gumbel2_pdf, ARGS3)
WRAP(gsl_cdf_gumbel2_P, ARGS3)
WRAP(gsl_cdf_gumbel2_Q, ARGS3)
WRAP(gsl_cdf_gumbel2_Pinv, ARGS3)
WRAP(gsl_cdf_gumbel2_Qinv, ARGS3)

#define WRAP_DISCRETE(func, args, uint_arg_names) \
  static double ampl##func(arglist *al) { \
    const void *has_names = uint_arg_names; \
    if (!check_args(al) || !check_uint_arg(al, 0, "k")) \
      return 0; \
    if (has_names) { \
      int i = 1; \
      for ( ; i < al->n; ++i) { \
        if (uint_arg_names[i] && !check_uint_arg(al, i, uint_arg_names[i])) \
          return 0; \
      } \
    } \
    if (al->derivs) \
      deriv_error(al, DERIVS_NOT_PROVIDED); \
    return check_result(al, func((unsigned)args)); \
  }

const char *const *const DEFAULT_ARGS = 0;

WRAP(gsl_ran_poisson, RNG_ARGS1)
WRAP_DISCRETE(gsl_ran_poisson_pdf, ARGS2, DEFAULT_ARGS)
WRAP_DISCRETE(gsl_cdf_poisson_P, ARGS2, DEFAULT_ARGS)
WRAP_DISCRETE(gsl_cdf_poisson_Q, ARGS2, DEFAULT_ARGS)

WRAP(gsl_ran_bernoulli, RNG_ARGS1)
WRAP_DISCRETE(gsl_ran_bernoulli_pdf, ARGS2, DEFAULT_ARGS)

static double amplgsl_ran_binomial(arglist *al) {
  if (!check_args(al) || !check_uint_arg(al, 1, "n"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al,
      gsl_ran_binomial(rng, al->ra[0], (unsigned)al->ra[1]));
}

const char *const BINOMIAL_ARGNAMES[] = {0, 0, "n"};

#define BINOMIAL_ARGS ARGS2, (unsigned)al->ra[2]

WRAP_DISCRETE(gsl_ran_binomial_pdf, BINOMIAL_ARGS, BINOMIAL_ARGNAMES)
WRAP_DISCRETE(gsl_cdf_binomial_P, BINOMIAL_ARGS, BINOMIAL_ARGNAMES)
WRAP_DISCRETE(gsl_cdf_binomial_Q, BINOMIAL_ARGS, BINOMIAL_ARGNAMES)

WRAP(gsl_ran_negative_binomial, RNG_ARGS2)
WRAP_DISCRETE(gsl_ran_negative_binomial_pdf, ARGS3, DEFAULT_ARGS)
WRAP_DISCRETE(gsl_cdf_negative_binomial_P, ARGS3, DEFAULT_ARGS)
WRAP_DISCRETE(gsl_cdf_negative_binomial_Q, ARGS3, DEFAULT_ARGS)

static double amplgsl_ran_pascal(arglist *al) {
  if (!check_args(al) || !check_uint_arg(al, 1, "n"))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al, gsl_ran_pascal(rng, al->ra[0], (unsigned)al->ra[1]));
}

WRAP_DISCRETE(gsl_ran_pascal_pdf, BINOMIAL_ARGS, BINOMIAL_ARGNAMES)
WRAP_DISCRETE(gsl_cdf_pascal_P, BINOMIAL_ARGS, BINOMIAL_ARGNAMES)
WRAP_DISCRETE(gsl_cdf_pascal_Q, BINOMIAL_ARGS, BINOMIAL_ARGNAMES)

WRAP(gsl_ran_geometric, RNG_ARGS1)
WRAP_DISCRETE(gsl_ran_geometric_pdf, ARGS2, DEFAULT_ARGS)
WRAP_DISCRETE(gsl_cdf_geometric_P, ARGS2, DEFAULT_ARGS)
WRAP_DISCRETE(gsl_cdf_geometric_Q, ARGS2, DEFAULT_ARGS)

static double amplgsl_ran_hypergeometric(arglist *al) {
  if (!check_args(al) || !check_uint_arg(al, 0, "n1") ||
      !check_uint_arg(al, 1, "n2") || !check_uint_arg(al, 2, "t")) {
    return 0;
  }
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  return check_result(al, gsl_ran_hypergeometric(rng,
      (unsigned)al->ra[0], (unsigned)al->ra[1], (unsigned)al->ra[2]));
}

const char *const HYPERGEOMETRIC_ARGNAMES[] = {0, "n1", "n2", "t"};

#define HYPERGEOMETRIC_ARGS \
  al->ra[0], (unsigned)al->ra[1], (unsigned)al->ra[2], (unsigned)al->ra[3]

WRAP_DISCRETE(gsl_ran_hypergeometric_pdf,
    HYPERGEOMETRIC_ARGS, HYPERGEOMETRIC_ARGNAMES)
WRAP_DISCRETE(gsl_cdf_hypergeometric_P,
    HYPERGEOMETRIC_ARGS, HYPERGEOMETRIC_ARGNAMES)
WRAP_DISCRETE(gsl_cdf_hypergeometric_Q,
    HYPERGEOMETRIC_ARGS, HYPERGEOMETRIC_ARGNAMES)

WRAP(gsl_ran_logarithmic, RNG_ARGS1)
WRAP_DISCRETE(gsl_ran_logarithmic_pdf, ARGS2, DEFAULT_ARGS)

#define ADDFUNC(name, num_args) \
    addfunc(#name, ampl##name, FUNCADD_REAL_VALUED, num_args, \
            const_cast<char*>(#name));

#define ADDFUNC_RANDOM(name, num_args) \
    addfunc(#name, ampl##name, FUNCADD_RANDOM_VALUED, num_args, \
            const_cast<char*>(#name));


// Mean, Standard Deviation and Variance
WRAP(gsl_stats_mean, ARRAY_ARGS);
WRAP(gsl_stats_variance, ARRAY_ARGS);
WRAP(gsl_stats_variance_m, ARRAY_ARGS_AND_LAST);
WRAP(gsl_stats_sd, ARRAY_ARGS);
WRAP(gsl_stats_sd_m, ARRAY_ARGS_AND_LAST);
WRAP(gsl_stats_tss, ARRAY_ARGS);
WRAP(gsl_stats_tss_m, ARRAY_ARGS_AND_LAST);
WRAP(gsl_stats_variance_with_fixed_mean, ARRAY_ARGS_AND_LAST);
WRAP(gsl_stats_sd_with_fixed_mean, ARRAY_ARGS_AND_LAST);

// Absolute deviation
WRAP(gsl_stats_absdev, ARRAY_ARGS);
WRAP(gsl_stats_absdev_m, ARRAY_ARGS_AND_LAST);

// Higher moments (skewness and kurtosis)
WRAP(gsl_stats_skew, ARRAY_ARGS);
WRAP(gsl_stats_skew_m_sd, ARRAY_ARGS_AND_LASTTWO);
WRAP(gsl_stats_kurtosis, ARRAY_ARGS);
WRAP(gsl_stats_kurtosis_m_sd, ARRAY_ARGS_AND_LASTTWO);

// Autocorrelation
WRAP(gsl_stats_lag1_autocorrelation, ARRAY_ARGS);
WRAP(gsl_stats_lag1_autocorrelation_m, ARRAY_ARGS_AND_LAST);


// Covariance
static double amplgsl_stats_covariance(arglist* al) {
  if (!check_args(al)) 
    return 0; 
  if (al->derivs) 
      deriv_error(al, DERIVS_NOT_PROVIDED); 
  size_t n = al->n / 2;
  double* p1 = al->ra;
  double* p2 = &al->ra[n];
  return check_result(al, gsl_stats_covariance(p1, 1, p2, 1, n));
}
static double amplgsl_stats_covariance_m(arglist* al) {
  if (!check_args(al))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  size_t n = (al->n - 2) / 2;
  double* p1 = al->ra;
  double* p2 = &al->ra[n];
  double m1 = al->ra[n * 2];
  double m2 = al->ra[n * 2 + 1];
  return check_result(al, gsl_stats_covariance_m(p1, 1, p2, 1, n, m1, m2));
}

// Correlation
static double amplgsl_stats_correlation(arglist* al) {
  if (!check_args(al))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  size_t n = al->n / 2;
  double* p1 = al->ra;
  double* p2 = &al->ra[n];
  return check_result(al, gsl_stats_correlation(p1, 1, p2, 1, n));
}

static double amplgsl_stats_spearman(arglist* al) {
  if (!check_args(al))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  size_t n = al->n / 2;
  double* p1 = al->ra;
  double* p2 = &al->ra[n];
  double* work = allocate_double(al, n * 2);
  return check_result(al, gsl_stats_spearman(p1, 1, p2, 1, n, work));
}

// Maximum and minimum
WRAP(gsl_stats_max, ARRAY_ARGS);
WRAP(gsl_stats_min, ARRAY_ARGS);
WRAP(gsl_stats_max_index, ARRAY_ARGS);
WRAP(gsl_stats_min_index, ARRAY_ARGS);

// Median and Percentiles
WRAP(gsl_stats_median_from_sorted_data, ARRAY_ARGS);
WRAP(gsl_stats_median, ARRAY_ARGS);
WRAP(gsl_stats_quantile_from_sorted_data, ARRAY_ARGS_AND_LAST);

// Order Statistics
WRAP(gsl_stats_select, ARRAY_ARGS_AND_LAST);

// Robust Location Estimates
WRAP(gsl_stats_trmean_from_sorted_data, ARRAY_ARGS_AND_LAST_TO_FIRST);
WRAP(gsl_stats_gastwirth_from_sorted_data, ARRAY_ARGS);



typedef double (*funcWithWork)(const double[], size_t, size_t, double[]);
static double impl_callWithAllocation(arglist* al, funcWithWork f) {
  if (!check_args(al))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  double* work = allocate_double(al, al->n);
  return check_result(al, f(al->ra, 1, al->n, work));
}

static double amplgsl_stats_mad0(arglist* al) {
  funcWithWork f = &gsl_stats_mad0;
  return impl_callWithAllocation(al,f);
}
static double amplgsl_stats_mad(arglist* al) {
  funcWithWork f = &gsl_stats_mad;
  return impl_callWithAllocation(al, f);
}
static double amplgsl_stats_Sn0_from_sorted_data(arglist* al) {
  funcWithWork f = &gsl_stats_Sn0_from_sorted_data;
  return impl_callWithAllocation(al, f);
}
static double amplgsl_stats_Sn_from_sorted_data(arglist* al) {
  funcWithWork f = &gsl_stats_Sn_from_sorted_data;
  return impl_callWithAllocation(al, f);
}
typedef double (*funcWithTwoWorks)(const double[], size_t, size_t, double[], int[]);
static double impl_callWithTwoAllocations(arglist* al, funcWithTwoWorks f) {
  if (!check_args(al))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  double* work = allocate_double(al, al->n);
  int* work_int = new int[al->n];
  double result = f(al->ra, 1, al->n, work, work_int);
  delete[] work_int;
  return check_result(al, result);
}
static double amplgsl_stats_Qn0_from_sorted_data(arglist* al) {
  funcWithTwoWorks f = &gsl_stats_Qn0_from_sorted_data;
  return impl_callWithTwoAllocations(al, f);
}
static double amplgsl_stats_Qn_from_sorted_data(arglist* al) {
  funcWithTwoWorks f = &gsl_stats_Qn_from_sorted_data;
  return impl_callWithTwoAllocations(al, f);
}

#include "gsl/gsl_sort.h"
static double amplgsl_sort(arglist* al) {
  if (!check_args(al))
    return 0;
  if (al->derivs)
    deriv_error(al, DERIVS_NOT_PROVIDED);
  gsl_sort(al->ra, 1, al->n);
  return 0;
}

extern "C" void funcadd_ASL(AmplExports *ae) {
  /* Don't call abort on error. */
  gsl_set_error_handler_off();
  
  addfunc("gsl_version", (rfunc)amplgsl_version,
      FUNCADD_STRING_VALUED, 0, const_cast<char*>("gsl_version"));
  /**
   * @file sorting
   *
   * Sorting Functions
   * -----------------
  */

  /**
  * .. function:: gsl_sort(x)
  *
  * This function sorts the input parameters.
  * The function can be used in AMPL to sort a variable or a parameter,
  * but it has one important limitation: it has to be declared on the 
  * same indexing expression that defines the entity to be ordered.
  * For example::
  * 
  *   set A;
  *   param p{A};
  *   function gsl_sort({a in A}(INOUT));
  *   call gsl_sort({a in A} d[a]); 
  * 
  */
  addfunc("gsl_sort", amplgsl_sort,
    FUNCADD_OUTPUT_ARGS, -1, const_cast<char*>("gsl_sort"));

  /**
   * @file elementary
   *
   * Elementary Functions
   * --------------------
   */

  /**
   * .. function:: gsl_log1p(x)
   *
   *  This function computes the value of $\log(1+x)$ in a way that is
   *  accurate for small $x$. It provides an alternative to the BSD math
   *  function ``log1p(x)``.
   */
  ADDFUNC(gsl_log1p, 1);

  /**
   * .. function:: gsl_expm1(x)
   *
   *  This function computes the value of $\exp(x)-1$ in a way that is
   *  accurate for small $x$. It provides an alternative to the BSD math
   *  function ``expm1(x)``.
   */
  ADDFUNC(gsl_expm1, 1);

  /**
   * .. function:: gsl_hypot(x, y)
   *
   *  This function computes the value of $\sqrt{x^2 + y^2}$ in a way that
   *  avoids overflow. It provides an alternative to the BSD math function
   *  ``hypot(x,y)``.
   */
  ADDFUNC(gsl_hypot, 2);

  /**
   * .. function:: gsl_hypot3(x, y, z)
   *
   *  This function computes the value of $\sqrt{x^2 + y^2 + z^2}$ in a way
   *  that avoids overflow.
   */
  ADDFUNC(gsl_hypot3, 3);

  /**
   * @file special
   *
   * Special Functions
   * -----------------
   *
   * This chapter describes the AMPL bindings for the GSL special function
   * library. The library includes routines for calculating the values of
   * Airy functions, Bessel functions, Clausen functions, Coulomb wave
   * functions, Coupling coefficients, the Dawson function, Debye functions,
   * Dilogarithms, Elliptic integrals, Jacobi elliptic functions, Error
   * functions, Exponential integrals, Fermi-Dirac functions, Gamma functions,
   * Gegenbauer functions, Hypergeometric functions, Laguerre functions,
   * Legendre functions and Spherical Harmonics, the Psi (Digamma) Function,
   * Synchrotron functions, Transport functions and Zeta functions.
   *
   * .. toctree::
   *    :maxdepth: 2
   *
   *    airy
   *    bessel
   *    clausen
   *    coulomb
   *    coupling
   *    dawson
   *    debye
   *    dilog
   *    ellint
   *    erf
   *    expint
   *    fermi-dirac
   *    gamma-beta
   *    gegenpoly
   *    hyperg
   *    laguerre
   *    lambert
   *    legendre
   *    log
   *    mathieu
   *    pow
   *    psi
   *    synchrotron
   *    transport
   *    zeta
   *    sf-refs
   */

  /* AMPL has built-in functions acosh, asinh and atanh so wrappers
     are not provided for their GSL equivalents. */

  /* Wrappers for functions operating on complex numbers are not provided
     since this requires support for structures/tuples as function arguments. */

  /**
   * @file airy
   *
   * Airy Functions and Derivatives
   * ==============================
   *
   * The Airy functions $\operatorname{Ai}(x)$ and $\operatorname{Bi}(x)$ are
   * defined by the integral representations,
   *
   * .. math::
   *   \operatorname{Ai}(x) = \frac{1}{\pi} \int_0^\infty
   *     \cos(\frac{1}{3} t^3 + xt) dt \\
   *   \operatorname{Bi}(x) = \frac{1}{\pi} \int_0^\infty
   *     (e^{-\frac{1}{3} t^3 + xt} + \sin(\frac{1}{3} t^3 + xt)) dt
   *
   * .. index:: Airy function
   *
   * For further information see Abramowitz & Stegun, Section 10.4.
   */

  /**
   * Airy Functions
   * --------------
   */

  /**
   * .. function:: gsl_sf_airy_Ai(x)
   *
   *  This routine computes the Airy function $\operatorname{Ai}(x)$.
   */
  ADDFUNC(gsl_sf_airy_Ai, 1);

  /**
   * .. function:: gsl_sf_airy_Bi(x)
   *
   *  This routine computes the Airy function $\operatorname{Bi}(x)$.
   */
  ADDFUNC(gsl_sf_airy_Bi, 1);

  /**
   * .. function:: gsl_sf_airy_Ai_scaled(x)
   *
   *  This routine computes a scaled version of the Airy function
   *  $\operatorname{S_A}(x) \operatorname{Ai}(x)$. For $x > 0$ the scaling
   *  factor $\operatorname{S_A}(x)$ is $\exp(+(2/3) x^{3/2})$, and is $1$
   *  for $x < 0$.
   */
  ADDFUNC(gsl_sf_airy_Ai_scaled, 1);

  /**
   * .. function:: gsl_sf_airy_Bi_scaled(x)
   *
   *  This routine computes a scaled version of the Airy function
   *  $\operatorname{S_B}(x) \operatorname{Bi}(x)$. For $x > 0$ the scaling
   *  factor $\operatorname{S_B}(x)$ is $\exp(-(2/3) x^{3/2})$, and is $1$
   *  for $x < 0$.
   */
  ADDFUNC(gsl_sf_airy_Bi_scaled, 1);

  /**
   * Zeros of Airy Functions
   * -----------------------
   */

  /**
   * .. function:: gsl_sf_airy_zero_Ai(s)
   *
   *  This routine computes the location of the $s$-th zero of the Airy
   *  function $\operatorname{Ai}(x)$.
   */
  ADDFUNC(gsl_sf_airy_zero_Ai, 1);

  /**
   * .. function:: gsl_sf_airy_zero_Bi(s)
   *
   *  This routine computes the location of the $s$-th zero of the Airy
   *  function $\operatorname{Bi}(x)$.
   */
  ADDFUNC(gsl_sf_airy_zero_Bi, 1);

  /**
   * Zeros of Derivatives of Airy Functions
   * --------------------------------------
   */

  /**
   * .. function:: gsl_sf_airy_zero_Ai_deriv(s)
   *
   *  This routine computes the location of the $s$-th zero of the Airy
   *  function derivative $\operatorname{Ai}'(x)$.
   */
  ADDFUNC(gsl_sf_airy_zero_Ai_deriv, 1);

  /**
   * .. function:: gsl_sf_airy_zero_Bi_deriv(s)
   *
   *  This routine computes the location of the $s$-th zero of the Airy
   *  function derivative $\operatorname{Bi}'(x)$.
   */
  ADDFUNC(gsl_sf_airy_zero_Bi_deriv, 1);

  /**
   * @file bessel
   *
   * Bessel Functions
   * ================
   *
   * The routines described in this section compute the Cylindrical Bessel
   * functions $J_n(x)$, $Y_n(x)$, Modified cylindrical Bessel functions
   * $I_n(x)$, $K_n(x)$, Spherical Bessel functions $j_l(x)$, $y_l(x)$,
   * and Modified Spherical Bessel functions $i_l(x)$, $k_l(x)$.
   * For more information see Abramowitz & Stegun, Chapters 9 and 10.
   *
   * .. index:: Bessel function
   */

  /**
   * Regular Cylindrical Bessel Functions
   * ------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_J0(x)
   *
   *  This routine computes the regular cylindrical Bessel function of
   *  zeroth order, $J_0(x)$.
   */
  ADDFUNC(gsl_sf_bessel_J0, 1);

  /**
   * .. function:: gsl_sf_bessel_J1(x)
   *
   *  This routine computes the regular cylindrical Bessel function of
   *  first order, $J_1(x)$.
   */
  ADDFUNC(gsl_sf_bessel_J1, 1);

  /**
   * .. function:: gsl_sf_bessel_Jn(n, x)
   *
   *  This routine computes the regular cylindrical Bessel function of
   *  integer order $n$, $J_n(x)$.
   */
  ADDFUNC(gsl_sf_bessel_Jn, 2);

  /**
   * Irregular Cylindrical Bessel Functions
   * --------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_Y0(x)
   *
   *  This routine computes the irregular cylindrical Bessel function of
   *  zeroth order, $Y_0(x)$, for $x >0$.
   */
  ADDFUNC(gsl_sf_bessel_Y0, 1);

  /**
   * .. function:: gsl_sf_bessel_Y1(x)
   *
   *  This routine computes the irregular cylindrical Bessel function of
   *  first order, $Y_1(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_Y1, 1);

  /**
   * .. function:: gsl_sf_bessel_Yn(n, x)
   *
   *  This routine computes the irregular cylindrical Bessel function of
   *  integer order $n$, $Y_n(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_Yn, 2);

  /**
   * Regular Modified Cylindrical Bessel Functions
   * ---------------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_I0(x)
   *
   *  This routine computes the regular modified cylindrical Bessel function
   *  of zeroth order, $I_0(x)$.
   */
  ADDFUNC(gsl_sf_bessel_I0, 1);

  /**
   * .. function:: gsl_sf_bessel_I1(x)
   *
   *  This routine computes the regular modified cylindrical Bessel function
   *  of first order, $I_1(x)$.
   */
  ADDFUNC(gsl_sf_bessel_I1, 1);

  /**
   * .. function:: gsl_sf_bessel_In(n, x)
   *
   *  This routine computes the regular modified cylindrical Bessel function
   *  of integer order $n$, $I_n(x)$.
   */
  ADDFUNC(gsl_sf_bessel_In, 2);

  /**
   * .. function:: gsl_sf_bessel_I0_scaled(x)
   *
   *  This routine computes the scaled regular modified cylindrical
   *  Bessel function of zeroth order $\exp(-|x|) I_0(x)$.
   */
  ADDFUNC(gsl_sf_bessel_I0_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_I1_scaled(x)
   *
   *  This routine computes the scaled regular modified cylindrical
   *  Bessel function of first order $\exp(-|x|) I_1(x)$.
   */
  ADDFUNC(gsl_sf_bessel_I1_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_In_scaled(n, x)
   *
   *  This routine computes the scaled regular modified cylindrical
   *  Bessel function of integer order $n$, $\exp(-|x|) I_n(x)$.
   */
  ADDFUNC(gsl_sf_bessel_In_scaled, 2);

  /**
   * Irregular Modified Cylindrical Bessel Functions
   * -----------------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_K0(x)
   *
   *  This routine computes the irregular modified cylindrical Bessel
   *  function of zeroth order, $K_0(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_K0, 1);

  /**
   * .. function:: gsl_sf_bessel_K1(x)
   *
   *  This routine computes the irregular modified cylindrical Bessel
   *  function of first order, $K_1(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_K1, 1);

  /**
   * .. function:: gsl_sf_bessel_Kn(n, x)
   *
   *  This routine computes the irregular modified cylindrical Bessel
   *  function of integer order $n$, $K_n(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_Kn, 2);

  /**
   * .. function:: gsl_sf_bessel_K0_scaled(x)
   *
   *  This routine computes the scaled irregular modified cylindrical Bessel
   *  function of zeroth order, $\exp(x) K_0(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_K0_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_K1_scaled(x)
   *
   *  This routine computes the scaled irregular modified cylindrical Bessel
   *  function of first order, $\exp(x) K_1(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_K1_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_Kn_scaled(n, x)
   *
   *  This routine computes the scaled irregular modified cylindrical Bessel
   *  function of integer order $n$, $\exp(x) K_n(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_Kn_scaled, 2);

  /**
   * Regular Spherical Bessel Functions
   * ----------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_j0(x)
   *
   *  This routine computes the regular spherical Bessel function of zeroth
   *  order, $j_0(x) = \sin(x)/x$.
   */
  ADDFUNC(gsl_sf_bessel_j0, 1);

  /**
   * .. function:: gsl_sf_bessel_j1(x)
   *
   *  This routine computes the regular spherical Bessel function of first
   *  order, $j_1(x) = (\sin(x)/x - \cos(x))/x$.
   */
  ADDFUNC(gsl_sf_bessel_j1, 1);

  /**
   * .. function:: gsl_sf_bessel_j2(x)
   *
   *  This routine computes the regular spherical Bessel function of second
   *  order, $j_2(x) = ((3/x^2 - 1)\sin(x) - 3\cos(x)/x)/x$.
   */
  ADDFUNC(gsl_sf_bessel_j2, 1);

  /**
   * .. function:: gsl_sf_bessel_jl(l, x)
   *
   *  This routine computes the regular spherical Bessel function of integer
   *  order $l$, $j_l(x)$, for $l \geq 0$ and $x \geq 0$.
   */
  ADDFUNC(gsl_sf_bessel_jl, 2);

  /**
   * Irregular Spherical Bessel Functions
   * ------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_y0(x)
   *
   *  This routine computes the irregular spherical Bessel function of
   *  zeroth order, $y_0(x) = -\cos(x)/x$.
   */
  ADDFUNC(gsl_sf_bessel_y0, 1);

  /**
   * .. function:: gsl_sf_bessel_y1(x)
   *
   *  This routine computes the irregular spherical Bessel function of
   *  first order, $y_1(x) = -(\cos(x)/x + \sin(x))/x$.
   */
  ADDFUNC(gsl_sf_bessel_y1, 1);

  /**
   * .. function:: gsl_sf_bessel_y2(x)
   *
   *  This routine computes the irregular spherical Bessel function of
   *  second order, $y_2(x) = (-3/x^3 + 1/x)\cos(x) - (3/x^2)\sin(x)$.
   */
  ADDFUNC(gsl_sf_bessel_y2, 1);

  /**
   * .. function:: gsl_sf_bessel_yl(l, x)
   *
   *  This routine computes the irregular spherical Bessel function of
   *  integer order $l$, $y_l(x)$, for $l \geq 0$.
   */
  ADDFUNC(gsl_sf_bessel_yl, 2);

  /**
   * Regular Modified Spherical Bessel Functions
   * -------------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_i0_scaled(x)
   *
   *  This routine computes the regular modified spherical Bessel function
   *  of zeroth order, $\exp(-|x|) i_0(x)$.
   */
  ADDFUNC(gsl_sf_bessel_i0_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_i1_scaled(x)
   *
   *  This routine computes the regular modified spherical Bessel function
   *  of first order, $\exp(-|x|) i_1(x)$.
   */
  ADDFUNC(gsl_sf_bessel_i1_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_i2_scaled(x)
   *
   *  This routine computes the regular modified spherical Bessel function
   *  of second order, $\exp(-|x|) i_2(x)$.
   */
  ADDFUNC(gsl_sf_bessel_i2_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_il_scaled(l, x)
   *
   *  This routine computes the regular modified spherical Bessel function
   *  of integer order $l$, $\exp(-|x|) i_l(x)$.
   */
  ADDFUNC(gsl_sf_bessel_il_scaled, 2);

  /**
   * Irregular Modified Spherical Bessel Functions
   * ---------------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_k0_scaled(x)
   *
   *  This routine computes the scaled irregular modified spherical Bessel
   *  function of zeroth order, $\exp(x) k_0(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_k0_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_k1_scaled(x)
   *
   *  This routine computes the scaled irregular modified spherical Bessel
   *  function of first order, $\exp(x) k_1(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_k1_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_k2_scaled(x)
   *
   *  This routine computes the scaled irregular modified spherical Bessel
   *  function of second order, $\exp(x) k_2(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_k2_scaled, 1);

  /**
   * .. function:: gsl_sf_bessel_kl_scaled(l, x)
   *
   *  This routine computes the scaled irregular modified spherical Bessel
   *  function of integer order $l$, $\exp(x) k_l(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_bessel_kl_scaled, 2);

  /**
   * Regular Bessel Function - Fractional Order
   * ------------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_Jnu(nu, x)
   *
   *  This routine computes the regular cylindrical Bessel function of
   *  fractional order $\nu$, $J_\nu(x)$.
   */
  ADDFUNC(gsl_sf_bessel_Jnu, 2);

  /**
   * Irregular Bessel Function - Fractional Order
   * ---------------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_Ynu(nu, x)
   *
   *  This routine computes the irregular cylindrical Bessel function of
   *  fractional order $\nu$, $Y_\nu(x)$.
   */
  ADDFUNC(gsl_sf_bessel_Ynu, 2);

  /**
   * Regular Modified Bessel Functions - Fractional Order
   * ----------------------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_Inu(nu, x)
   *
   *  This routine computes the regular modified Bessel function of
   *  fractional order $\nu$, $I_\nu(x)$ for $x > 0$, $\nu > 0$.
   */
  ADDFUNC(gsl_sf_bessel_Inu, 2);

  /**
   * .. function:: gsl_sf_bessel_Inu_scaled(nu, x)
   *
   *  This routine computes the scaled regular modified Bessel function of
   *  fractional order $\nu$, $\exp(-|x|) I_\nu(x)$ for $x > 0$,
   *  $\nu > 0$.
   */
  ADDFUNC(gsl_sf_bessel_Inu_scaled, 2);

  /**
   * Irregular Modified Bessel Functions - Fractional Order
   * ------------------------------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_Knu(nu, x)
   *
   *  This routine computes the irregular modified Bessel function of
   *  fractional order $\nu$, $K_\nu(x)$ for $x > 0$, $\nu > 0$.
   */
  ADDFUNC(gsl_sf_bessel_Knu, 2);

  /**
   * .. function:: gsl_sf_bessel_lnKnu(nu, x)
   *
   *  This routine computes the logarithm of the irregular modified Bessel
   *  function of fractional order $\nu$, $\ln(K_\nu(x))$ for $x > 0$,
   *  $\nu > 0$.
   */
  ADDFUNC(gsl_sf_bessel_lnKnu, 2);

  /**
   * .. function:: gsl_sf_bessel_Knu_scaled(nu, x)
   *
   *  This routine computes the scaled irregular modified Bessel function of
   *  fractional order $\nu$, $\exp(|x|) K_\nu(x)$ for $x > 0$, $\nu > 0$.
   */
  ADDFUNC(gsl_sf_bessel_Knu_scaled, 2);

  /**
   * Zeros of Regular Bessel Functions
   * ---------------------------------
   */

  /**
   * .. function:: gsl_sf_bessel_zero_J0(s)
   *
   *  This routine computes the location of the $s$-th positive zero of the
   *  Bessel function $J_0(x)$.
   */
  ADDFUNC(gsl_sf_bessel_zero_J0, 1);

  /**
   * .. function:: gsl_sf_bessel_zero_J1(s)
   *
   *  This routine computes the location of the $s$-th positive zero of the
   *  Bessel function $J_1(x)$.
   */
  ADDFUNC(gsl_sf_bessel_zero_J1, 1);

  /**
   * .. function:: gsl_sf_bessel_zero_Jnu(nu, s)
   *
   *  This routine computes the location of the $s$-th positive zero of the
   *  Bessel function $J_\nu(x)$. The current implementation does not support
   *  negative values of ``nu``.
   */
  ADDFUNC(gsl_sf_bessel_zero_Jnu, 2);

  /**
   * @file clausen
   *
   * Clausen Function
   * ================
   *
   * The Clausen function is defined by the following integral,
   *
   * .. math::
   *  \operatorname{Cl_2}(x) = -\int_0^x \log(2 \sin(t/2)) dt
   *
   * .. index:: Clausen function
   *
   * It is related to the dilogarithm by
   *
   * .. math::
   *   \operatorname{Cl_2}(\theta) =
   *     \operatorname{Im} \operatorname{Li_2}(\exp(i\theta)).
   */

  /**
   * .. function:: gsl_sf_clausen(x)
   *
   *  This routine computes the :index:`Clausen integral`
   *  $\operatorname{Cl_2}(x)$.
   */
  ADDFUNC(gsl_sf_clausen, 1);

  /**
   * @file coulomb
   *
   * Coulomb Functions
   * =================
   *
   * .. index:: Coulomb function
   */

  /**
   * Normalized Hydrogenic Bound States
   * ----------------------------------
   */

  /**
   * .. function:: gsl_sf_hydrogenicR_1(Z, r)
   *
   *  This routine computes the lowest-order normalized hydrogenic bound
   *  state radial wavefunction $R_1 := 2 Z \sqrt{Z} \exp(-Z r)$.
   */
  ADDFUNC(gsl_sf_hydrogenicR_1, 2);

  /**
   * .. function:: gsl_sf_hydrogenicR(n, l, Z, r)
   *
   *  This routine computes the $n$-th normalized hydrogenic bound state
   *  radial wavefunction,
   *
   *  .. math::
   *    R_n := 2 (Z^{3/2}/n^2) \sqrt{(n-l-1)!/(n+l)!} \exp(-Z r/n) (2Zr/n)^l
   *              L^{2l+1}_{n-l-1}(2Zr/n).
   *
   *  where $L^a_b(x)$ is the generalized Laguerre polynomial
   *  (see :ref:`laguerre-functions`). The normalization is chosen such that
   *  the wavefunction $\psi$ is given by $\psi(n,l,r) = R_n Y_{lm}$.
   */
  ADDFUNC(gsl_sf_hydrogenicR, 4);

  /**
   * Coulomb Wave Function Normalization Constant
   * --------------------------------------------
   *
   * The Coulomb wave function normalization constant is defined in
   * Abramowitz 14.1.7.
   */

  /**
   * .. function:: gsl_sf_coulomb_CL(L, eta)
   *
   *  This function computes the Coulomb wave function normalization
   *  constant $C_L(\eta)$ for $L > -1$.
   */
  ADDFUNC(gsl_sf_coulomb_CL, 2);

  /**
   * @file coupling
   *
   * Coupling Coefficients
   * =====================
   *
   * The Wigner 3-j, 6-j and 9-j symbols give the coupling coefficients
   * for combined angular momentum vectors. Since the arguments of the
   * standard coupling coefficient functions are integer or half-integer,
   * the arguments of the following functions are, by convention,
   * integers equal to twice the actual spin value. For information on
   * the 3-j coefficients see Abramowitz & Stegun, Section 27.9.
   *
   * .. index:: coupling coefficient
   */

  /**
   * .. function:: gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc)
   *
   *  These routines compute the Wigner 3-j coefficient,
   *
   *  .. math::
   *    \left( \begin{array}{ccc}
   *           ja & jb & jc \\
   *           ma & mb & mc
   *           \end{array} \right)
   *
   *  where the arguments are given in half-integer units, ja = two_ja / 2,
   *  ma = two_ma / 2, etc.
   *
   * .. index:: Wigner 3-j coefficient
   */
  ADDFUNC(gsl_sf_coupling_3j, 6);

  /**
   * .. function::
   *   gsl_sf_coupling_6j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf)
   *
   *  These routines compute the Wigner 6-j coefficient,
   *
   *  .. math::
   *    \left( \begin{array}{ccc}
   *           ja & jb & jc \\
   *           jd & je & jf
   *           \end{array} \right)
   *
   *  where the arguments are given in half-integer units, ja = two_ja / 2,
   *  jb = two_jb / 2, etc.
   *
   * .. index:: Wigner 6-j coefficient
   */
  ADDFUNC(gsl_sf_coupling_6j, 6);

  /**
   * .. function::
   *   gsl_sf_coupling_9j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji)
   *
   *  These routines compute the Wigner 9-j coefficient,
   *
   *  .. math::
   *    \left( \begin{array}{ccc}
   *           ja & jb & jc \\
   *           jd & je & jf \\
   *           jg & jh & ji
   *           \end{array} \right)
   *
   *  where the arguments are given in half-integer units, ja = two_ja / 2,
   *  jb = two_jb / 2, etc.
   *
   * .. index:: Wigner 9-j coefficient
   */
  ADDFUNC(gsl_sf_coupling_9j, 9);

  /**
   * @file dawson
   *
   * Dawson Function
   * ===============
   *
   * .. index:: Dawson function
   *
   * The :index:`Dawson integral` is defined by
   * $\exp(-x^2) \int_0^x \exp(t^2) dt$.
   * A table of Dawson's integral can be found in Abramowitz & Stegun,
   * Table 7.5.
   */

  /**
   * .. function:: gsl_sf_dawson(x)
   *
   *  This routine computes the value of Dawson's integral for $x$.
   */
  ADDFUNC(gsl_sf_dawson, 1);

  /**
   * @file debye
   *
   * Debye Functions
   * ===============
   *
   * The Debye functions $D_n(x)$ are defined by the following integral,
   *
   * .. math::
   *   D_n(x) = n/x^n \int_0^x (t^n/(e^t - 1)) dt
   *
   * .. index:: Debye function
   *
   * For further information see Abramowitz & Stegun, Section 27.1.
   */

  /**
   * .. function:: gsl_sf_debye_1(x)
   *
   *  This routine computes the first-order Debye function
   *  $D_1(x) = (1/x) \int_0^x (t/(e^t - 1)) dt$.
   */
  ADDFUNC(gsl_sf_debye_1, 1);

  /**
   * .. function:: gsl_sf_debye_2(x)
   *
   *  This routine computes the second-order Debye function
   *  $D_2(x) = (2/x^2) \int_0^x (t^2/(e^t - 1)) dt$.
   */
  ADDFUNC(gsl_sf_debye_2, 1);

  /**
   * .. function:: gsl_sf_debye_3(x)
   *
   *  This routine computes the third-order Debye function
   *  $D_3(x) = (3/x^3) \int_0^x (t^3/(e^t - 1)) dt$.
   */
  ADDFUNC(gsl_sf_debye_3, 1);

  /**
   * .. function:: gsl_sf_debye_4(x)
   *
   *  This routine computes the fourth-order Debye function
   *  $D_4(x) = (4/x^4) \int_0^x (t^4/(e^t - 1)) dt$.
   */
  ADDFUNC(gsl_sf_debye_4, 1);

  /**
   * .. function:: gsl_sf_debye_5(x)
   *
   *  This routine computes the fifth-order Debye function
   *  $D_5(x) = (5/x^5) \int_0^x (t^5/(e^t - 1)) dt$.
   */
  ADDFUNC(gsl_sf_debye_5, 1);

  /**
   * .. function:: gsl_sf_debye_6(x)
   *
   *  This routine computes the sixth-order Debye function
   *  $D_6(x) = (6/x^6) \int_0^x (t^6/(e^t - 1)) dt$.
   */
  ADDFUNC(gsl_sf_debye_6, 1);

  /**
   * @file dilog
   *
   * Dilogarithm
   * ===========
   *
   * .. index:: dilogarithm
   */

  /**
   * .. function:: gsl_sf_dilog(x)
   *
   *  This routine computes the dilogarithm for a real argument. In Lewin's
   *  notation this is $\operatorname{Li}_2(x)$, the real part of the
   *  dilogarithm of a real $x$. It is defined by the integral representation
   *  $\operatorname{Li}_2(x) = -\operatorname{Re}\int_0^x (\log(1-s) / s) ds$.
   *  Note that $\operatorname{Im}(\operatorname{Li}_2(x)) = 0$ for $x <= 1$,
   *  and $-\pi\log(x)$ for $x > 1$.
   *
   *  Note that Abramowitz & Stegun refer to the Spence integral
   *  $S(x)=\operatorname{Li}_2(1-x)$ as the dilogarithm rather than
   *  $\operatorname{Li}_2(x)$.
   */
  ADDFUNC(gsl_sf_dilog, 1);

  /**
   * @file ellint
   *
   * Elliptic Integrals
   * ==================
   *
   * Information about the elliptic integrals can be found in
   * Abramowitz & Stegun, Chapter 17.
   *
   * .. index:: elliptic integral
   *
   * Definition of Legendre Forms
   * ----------------------------
   *
   * The Legendre forms of elliptic integrals $F(\phi,k)$, $E(\phi,k)$ and
   * $\Pi(\phi,k,n)$ are defined by,
   *
   * .. math::
   *   F(\phi,k) = \int_0^\phi 1/\sqrt{1 - k^2 \sin^2(t)} dt \\
   *   E(\phi,k) = \int_0^\phi \sqrt{1 - k^2 \sin^2(t)} dt \\
   *   \Pi(\phi,k,n) =
   *     \int_0^\phi 1/((1 + n \sin^2(t))\sqrt{1 - k^2 \sin^2(t)}) dt
   *
   * The complete Legendre forms are denoted by $K(k) = F(\pi/2, k)$ and
   * $E(k) = E(\pi/2, k)$.
   *
   * The notation used here is based on Carlson, Numerische Mathematik 33
   * (1979) 1 and differs slightly from that used by Abramowitz & Stegun,
   * where the functions are given in terms of the parameter $m = k^2$ and
   * $n$ is replaced by $-n$.
   *
   * Definition of Carlson Forms
   * ---------------------------
   *
   * The Carlson symmetric forms of elliptical integrals $RC(x,y)$,
   * $RD(x,y,z)$, $RF(x,y,z)$ and $RJ(x,y,z,p)$ are defined by,
   *
   * .. math::
   *   RC(x,y) = 1/2 \int_0^\infty (t+x)^{-1/2} (t+y)^{-1} dt \\
   *   RD(x,y,z) = 3/2 \int_0^\infty
   *     (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-3/2} dt \\
   *   RF(x,y,z) = 1/2 \int_0^\infty
   *     (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-1/2} dt \\
   *   RJ(x,y,z,p) = 3/2 \int_0^\infty
   *     (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-1/2} (t+p)^{-1} dt
   */

  /**
   * Legendre Form of Complete Elliptic Integrals
   * --------------------------------------------
   */

  /**
   * .. function:: gsl_sf_ellint_Kcomp(k)
   *
   *  This routine computes the complete elliptic integral $K(k)$.
   *  Note that Abramowitz & Stegun define this function in terms
   *  of the parameter $m = k^2$.
   */
  ADDFUNC(gsl_sf_ellint_Kcomp, 1);

  /**
   * .. function:: gsl_sf_ellint_Ecomp(k)
   *
   *  This routine computes the complete elliptic integral $E(k)$.
   *  Note that Abramowitz & Stegun define this function in terms
   *  of the parameter $m = k^2$.
   */
  ADDFUNC(gsl_sf_ellint_Ecomp, 1);

  /**
   * .. function:: gsl_sf_ellint_Pcomp(k, n)
   *
   *  This routine computes the complete elliptic integral $\Pi(k,n)$.
   *  Note that Abramowitz & Stegun define this function in terms
   *  of the parameters $m = k^2$ and $\sin^2(\alpha) = k^2$, with the
   *  change of sign $n \to -n$.
   */
  ADDFUNC(gsl_sf_ellint_Pcomp, 2);

  /**
   * Legendre Form of Incomplete Elliptic Integrals
   * ----------------------------------------------
   */

  /**
   * .. function:: gsl_sf_ellint_F(phi, k)
   *
   *  This routine computes the incomplete elliptic integral $F(\phi,k)$.
   *  Note that Abramowitz & Stegun define this function in terms of the
   *  parameter $m = k^2$.
   */
  ADDFUNC(gsl_sf_ellint_F, 2);

  /**
   * .. function:: gsl_sf_ellint_E(phi, k)
   *
   *  This routine computes the incomplete elliptic integral $E(\phi,k)$.
   *  Note that Abramowitz & Stegun define this function in terms of the
   *  parameter $m = k^2$.
   */
  ADDFUNC(gsl_sf_ellint_E, 2);

  /**
   * .. function:: gsl_sf_ellint_P(phi, k, n)
   *
   *  This routine computes the incomplete elliptic integral $\Pi(\phi,k,n)$.
   *  Note that Abramowitz & Stegun define this function in terms of the
   *  parameters $m = k^2$ and $\sin^2(\alpha) = k^2$, with the change of
   *  sign $n \to -n$.
   */
  ADDFUNC(gsl_sf_ellint_P, 3);

  /**
   * .. function:: gsl_sf_ellint_D(phi, k)
   *
   *  This routine computes the incomplete elliptic integral $D(\phi,k)$
   *  which is defined through the Carlson form $RD(x,y,z)$ by the following
   *  relation,
   *
   *  .. math::
   *    D(\phi,k) = (1/3)(\sin(\phi))^3
   *      RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1).
   */
  ADDFUNC(gsl_sf_ellint_D, 2);

  /**
   * Carlson Forms
   * -------------
   */

  /**
   * .. function:: gsl_sf_ellint_RC(x, y)
   *
   *  This routine computes the incomplete elliptic integral $RC(x,y)$.
   */
  ADDFUNC(gsl_sf_ellint_RC, 2);

  /**
   * .. function:: gsl_sf_ellint_RD(x, y, z)
   *
   *  This routine computes the incomplete elliptic integral $RD(x,y,z)$.
   */
  ADDFUNC(gsl_sf_ellint_RD, 3);

  /**
   * .. function:: gsl_sf_ellint_RF(x, y, z)
   *
   *  This routine computes the incomplete elliptic integral $RF(x,y,z)$.
   */
  ADDFUNC(gsl_sf_ellint_RF, 3);

  /**
   * .. function:: gsl_sf_ellint_RJ(x, y, z, p)
   *
   *  This routine computes the incomplete elliptic integral $RJ(x,y,z,p)$.
   */
  ADDFUNC(gsl_sf_ellint_RJ, 4);

  /* Elliptic Functions (Jacobi) */
  /* Wrapper for gsl_sf_elljac_e is not provided since the latter produces
     multiple values (through output parameters). */

  /**
   * @file erf
   *
   * Error Functions
   * ===============
   *
   * The :index:`error function` is described in Abramowitz & Stegun,
   * Chapter 7.
   */

  /**
   * .. function:: gsl_sf_erf(x)
   *
   *  This routine computes the error function $\operatorname{erf}(x)$, where
   *
   *  .. math::
   *    \operatorname{erf}(x) = (2/\sqrt{\pi}) \int_0^x \exp(-t^2) dt.
   */
  ADDFUNC(gsl_sf_erf, 1);

  /**
   * .. function:: gsl_sf_erfc(x)
   *
   *  This routine computes the complementary error function
   *
   *  .. math::
   *    \operatorname{erfc}(x) = 1 - \operatorname{erf}(x) =
   *      (2/\sqrt{\pi}) \int_x^\infty \exp(-t^2) dt.
   */
  ADDFUNC(gsl_sf_erfc, 1);

  /**
   * .. function:: gsl_sf_log_erfc(x)
   *
   *  This routine computes the logarithm of the complementary error function
   *  $\log(\operatorname{erfc}(x))$.
   */
  ADDFUNC(gsl_sf_log_erfc, 1);

  /**
   * Probability functions
   * ---------------------
   *
   * The probability functions for the Normal or Gaussian distribution are
   * described in Abramowitz & Stegun, Section 26.2.
   */

  /**
   * .. function:: gsl_sf_erf_Z(x)
   *
   *  This routine computes the Gaussian probability density function
   *  $Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2)$.
   */
  ADDFUNC(gsl_sf_erf_Z, 1);

  /**
   * .. function:: gsl_sf_erf_Q(x)
   *
   *  This routine computes the upper tail of the Gaussian probability
   *  function $Q(x) = (1/\sqrt{2\pi}) \int_x^\infty \exp(-t^2/2) dt$.
   */
  ADDFUNC(gsl_sf_erf_Q, 1);

  /**
   * The hazard function for the normal distribution, also known as the
   * inverse Mills' ratio, is defined as,
   *
   * .. math::
   *   h(x) = Z(x)/Q(x) =
   *     \sqrt{2/\pi} \exp(-x^2 / 2) / \operatorname{erfc}(x/\sqrt 2)
   *
   * It decreases rapidly as $x$ approaches $-\infty$ and asymptotes
   * to $h(x) \sim x$ as $x$ approaches $+\infty$.
   */

  /**
   * .. function:: gsl_sf_hazard(x)
   *
   *  This routine computes the hazard function for the normal distribution.
   */
  ADDFUNC(gsl_sf_hazard, 1);

  /**
   * @file expint
   *
   * Exponential Integrals
   * =====================
   *
   * .. index:: exponential integral
   */

  /**
   * Exponential Integral
   * --------------------
   */

  /**
   * .. function:: gsl_sf_expint_E1(x)
   *
   *  This routine computes the exponential integral $\operatorname{E_1}(x)$,
   *
   *  .. math::
   *    \operatorname{E_1}(x) :=
   *      \operatorname{Re} \int_1^\infty \exp(-xt)/t dt.
   */
  ADDFUNC(gsl_sf_expint_E1, 1);

  /**
   * .. function:: gsl_sf_expint_E2(x)
   *
   *  This routine computes the second-order exponential integral
   *  $\operatorname{E_2}(x)$,
   *
   *  .. math::
   *    \operatorname{E_2(x)} :=
   *      \operatorname{Re} \int_1^\infty \exp(-xt)/t^2 dt.
   */
  ADDFUNC(gsl_sf_expint_E2, 1);

  /**
   * .. function:: gsl_sf_expint_En(n, x)
   *
   *  This routine computes the exponential integral $\operatorname{E_n}(x)$
   *  of order $n$,
   *
   *  .. math::
   *    \operatorname{E_n}(x) :=
   *      \operatorname{Re} \int_1^\infty \exp(-xt)/t^n dt.
   */
  ADDFUNC(gsl_sf_expint_En, 2);

  /**
   * Ei(x)
   * -----
   */

  /**
   * .. function:: gsl_sf_expint_Ei(x)
   *
   *  These routines compute the exponential integral $\operatorname{Ei}(x)$,
   *
   *  .. math::
   *    \operatorname{Ei}(x) := - PV(\int_{-x}^\infty \exp(-t)/t dt)
   *
   *  where $PV$ denotes the principal value of the integral.
   */
  ADDFUNC(gsl_sf_expint_Ei, 1);

  /**
   * Hyperbolic Integrals
   * --------------------
   */

  /**
   * .. function:: gsl_sf_Shi(x)
   *
   *  This routine computes the integral
   *
   *  .. math::
   *    \operatorname{Shi}(x) = \int_0^x \sinh(t)/t dt.
   */
  ADDFUNC(gsl_sf_Shi, 1);

  /**
   * .. function:: gsl_sf_Chi(x)
   *
   *  This routine computes the integral
   *
   *  .. math::
   *    \operatorname{Chi}(x) := \operatorname{Re}[
   *      \gamma_E + \log(x) + \int_0^x (\cosh(t)-1)/t dt],
   *
   *  where $\gamma_E$ is the Euler constant.
   */
  ADDFUNC(gsl_sf_Chi, 1);

  /**
   * Ei_3(x)
   * -------
   */

  /**
   * .. function:: gsl_sf_expint_3(x)
   *
   *  This routine computes the third-order exponential integral
   *
   *  .. math::
   *    \operatorname{Ei_3}(x) = \int_0^x \exp(-t^3) dt \text{ for } x \geq 0.
   */
  ADDFUNC(gsl_sf_expint_3, 1);

  /**
   * Trigonometric Integrals
   * -----------------------
   */

  /**
   * .. function:: gsl_sf_Si(x)
   *
   *  This routine computes the :index:`Sine integral`
   *
   *  .. math::
   *    \operatorname{Si}(x) = \int_0^x \sin(t)/t dt.
   */
  ADDFUNC(gsl_sf_Si, 1);

  /**
   * .. function:: gsl_sf_Ci(x)
   *
   *  This routine computes the :index:`Cosine integral`
   *
   *  .. math::
   *    \operatorname{Ci}(x) = -\int_x^\infty \cos(t)/t dt \text{ for } x > 0.
   */
  ADDFUNC(gsl_sf_Ci, 1);

  /**
   * Arctangent Integral
   * -------------------
   */

  /**
   * .. function:: gsl_sf_atanint(x)
   *
   *  This routine computes the :index:`Arctangent integral`, which is
   *  defined as
   *
   *  .. math::
   *    \operatorname{AtanInt}(x) = \int_0^x \arctan(t)/t dt.
   */
  ADDFUNC(gsl_sf_atanint, 1);

  /**
   * @file fermi-dirac
   *
   * Fermi-Dirac Function
   * ====================
   *
   * .. index:: Fermi-Dirac function
   */

  /**
   * Complete Fermi-Dirac Integrals
   * ------------------------------
   *
   * The complete :index:`Fermi-Dirac integral` $F_j(x)$ is given by,
   *
   * .. math::
   *   F_j(x) := (1/\Gamma(j+1)) \int_0^\infty (t^j / (\exp(t-x) + 1)) dt
   *
   * Note that the Fermi-Dirac integral is sometimes defined without the
   * normalisation factor in other texts.
   */

  /**
   * .. function:: gsl_sf_fermi_dirac_m1(x)
   *
   *  This routine computes the complete Fermi-Dirac integral with an index
   *  of -1. This integral is given by $F_{-1}(x) = e^x / (1 + e^x)$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_m1, 1);

  /**
   * .. function:: gsl_sf_fermi_dirac_0(x)
   *
   *  This routine computes the complete Fermi-Dirac integral with an index
   *  of 0. This integral is given by $F_0(x) = \ln(1 + e^x)$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_0, 1);


  /**
   * .. function:: gsl_sf_fermi_dirac_1(x)
   *
   *  This routine computes the complete Fermi-Dirac integral with an index
   *  of 1, $F_1(x) = \int_0^\infty (t /(\exp(t-x)+1)) dt$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_1, 1);

  /**
   * .. function:: gsl_sf_fermi_dirac_2(x)
   *
   *  This routine computes the complete Fermi-Dirac integral with an index
   *  of 2, $F_2(x) = (1/2) \int_0^\infty (t^2 /(\exp(t-x)+1)) dt$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_2, 1);

  /**
   * .. function:: gsl_sf_fermi_dirac_int(j, x)
   *
   *  This routine computes the complete Fermi-Dirac integral with an
   *  integer index of $j$,
   *  $F_j(x) = (1/\Gamma(j+1)) \int_0^\infty (t^j /(\exp(t-x)+1)) dt$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_int, 2);

  /**
   * .. function:: gsl_sf_fermi_dirac_mhalf(x)
   *
   *  This routine computes the complete Fermi-Dirac integral $F_{-1/2}(x)$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_mhalf, 1);

  /**
   * .. function:: gsl_sf_fermi_dirac_half(x)
   *
   *  This routine computes the complete Fermi-Dirac integral $F_{1/2}(x)$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_half, 1);

  /**
   * .. function:: gsl_sf_fermi_dirac_3half(x)
   *
   *  This routine computes the complete Fermi-Dirac integral $F_{3/2}(x)$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_3half, 1);

  /**
   * Incomplete Fermi-Dirac Integrals
   * --------------------------------
   *
   * The incomplete Fermi-Dirac integral F_j(x,b) is given by,
   *
   * .. math::
   *   F_j(x,b) := (1/\Gamma(j+1)) \int_b^\infty (t^j / (\exp(t-x) + 1)) dt
   */

  /**
   * .. function:: gsl_sf_fermi_dirac_inc_0(x, b)
   *
   *  This routine computes the incomplete Fermi-Dirac integral with an index
   *  of zero, $F_0(x,b) = \ln(1 + e^{b-x}) - (b-x)$.
   */
  ADDFUNC(gsl_sf_fermi_dirac_inc_0, 2);

  /**
   * @file gamma-beta
   *
   * Gamma and Beta Functions
   * ========================
   *
   * This following routines compute the gamma and beta functions in their
   * full and incomplete forms.
   *
   * .. index:
   *   Gamma function
   *   Beta function
   */

  /**
   * .. _gamma-functions:
   *
   * Gamma Functions
   * ---------------
   *
   * The :index:`Gamma function` is defined by the following integral,
   *
   * .. math::
   *   \Gamma(x) = \int_0^\infty t^{x-1} \exp(-t) dt
   *
   * It is related to the factorial function by $\Gamma(n)=(n-1)!$ for
   * positive integer $n$. Further information on the Gamma function can
   * be found in Abramowitz & Stegun, Chapter 6.
   */

  /**
   * .. function:: gsl_sf_gamma(x)
   *
   *  This routine computes the Gamma function $\Gamma(x)$, subject to $x$
   *  not being a negative integer or zero. The function is computed using
   *  the real Lanczos method. The maximum value of $x$ such that $\Gamma(x)$
   *  is not considered an overflow is 171.0.
   */
  ADDFUNC(gsl_sf_gamma, 1);

  /**
   * .. function:: gsl_sf_lngamma(x)
   *
   *  This routine computes the logarithm of the Gamma function,
   *  $\log(\Gamma(x))$, subject to $x$ not being a negative integer or zero.
   *  For $x<0$ the real part of $\log(\Gamma(x))$ is returned, which is
   *  equivalent to $\log(|\Gamma(x)|)$. The function is computed using the
   *  real Lanczos method.
   */
  ADDFUNC(gsl_sf_lngamma, 1);

  /**
   * .. function:: gsl_sf_gammastar(x)
   *
   *  This routine computes the regulated Gamma Function $\Gamma^*(x)$ for
   *  $x > 0$. The regulated gamma function is given by,
   *
   *  .. math::
   *    \Gamma^*(x) = \Gamma(x)/(\sqrt{2\pi} x^{(x-1/2)} \exp(-x))
   *                = (1 + (1/12x) + ...)  \text{ for } x \to \infty
   *
   *  and is a useful suggestion of Temme.
   */
  ADDFUNC(gsl_sf_gammastar, 1);

  /**
   * .. function:: gsl_sf_gammainv(x)
   *
   *  This routine computes the reciprocal of the gamma function,
   *  $1/\Gamma(x)$ using the real Lanczos method.
   */
  ADDFUNC(gsl_sf_gammainv, 1);

  /* Wrapper for factorials are not provided since these are easily
     implemented using built-in AMPL features like the prod operator. */

  /**
   * Pochhammer Symbol
   * -----------------
   */

  /**
   * .. function:: gsl_sf_poch(a, x)
   *
   *  This routine computes the :index:`Pochhammer symbol`
   *  $(a)_x = \Gamma(a + x)/\Gamma(a)$. The Pochhammer symbol is also
   *  known as the Apell symbol and sometimes written as $(a,x)$.
   *  When $a$ and $a+x$ are negative integers or zero, the limiting
   *  value of the ratio is returned.
   */
  ADDFUNC(gsl_sf_poch, 2);

  /**
   * .. function:: gsl_sf_lnpoch(a, x)
   *
   *  This routine computes the logarithm of the Pochhammer symbol,
   *  $\log((a)_x) = \log(\Gamma(a + x)/\Gamma(a))$.
   */
  ADDFUNC(gsl_sf_lnpoch, 2);

  /**
   * .. function:: gsl_sf_pochrel(a, x)
   *
   *  This routine computes the relative Pochhammer symbol
   *  $((a)_x - 1)/x$ where $(a)_x = \Gamma(a + x)/\Gamma(a)$.
   */
  ADDFUNC(gsl_sf_pochrel, 2);

  /**
   * Incomplete Gamma Functions
   * --------------------------
   */

  /**
   * .. function:: gsl_sf_gamma_inc(a, x)
   *
   *  This routine computes the unnormalized incomplete :index:`Gamma Function`
   *  $\Gamma(a,x) = \int_x^\infty t^{a-1} \exp(-t) dt$ for a real and
   *  $x \geq 0$.
   */
  ADDFUNC(gsl_sf_gamma_inc, 2);

  /**
   * .. function:: gsl_sf_gamma_inc_Q(a, x)
   *
   *  This routine computes the normalized incomplete Gamma Function
   *  $Q(a,x) = 1/\Gamma(a) \int_x^\infty t^{a-1} \exp(-t) dt$ for
   *  $a > 0$, $x \geq 0$.
   */
  ADDFUNC(gsl_sf_gamma_inc_Q, 2);

  /**
   * .. function:: gsl_sf_gamma_inc_P(a, x)
   *
   *  This routine computes the complementary normalized incomplete
   *  Gamma Function
   *
   *  .. math::
   *    P(a,x) = 1 - Q(a,x) = 1/\Gamma(a) \int_0^x t^{a-1} \exp(-t) dt
   *    \text{ for } a > 0, x \geq 0.
   *
   *  Note that Abramowitz & Stegun call $P(a,x)$ the incomplete gamma
   *  function (section 6.5).
   */
  ADDFUNC(gsl_sf_gamma_inc_P, 2);

  /**
   * Beta Functions
   * --------------
   */

  /**
   * .. function:: gsl_sf_beta(a, b)
   *
   *  This routine computes the :index:`Beta Function`,
   *  $\operatorname{B}(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b)$
   *  subject to $a$ and $b$ not being negative integers.
   */
  ADDFUNC(gsl_sf_beta, 2);

  /**
   * .. function:: gsl_sf_lnbeta(a, b)
   *
   *  This routine computes the logarithm of the Beta Function,
   *  $\log(\operatorname{B}(a,b))$ subject to $a$ and $b$ not being
   *  negative integers.
   */
  ADDFUNC(gsl_sf_lnbeta, 2);

  /**
   * Incomplete Beta Function
   * ------------------------
   */

  /**
   * .. function:: gsl_sf_beta_inc(a, b, x)
   *
   *  This routine computes the normalized incomplete Beta function
   *  $I_x(a,b) = \operatorname{B}_x(a,b)/\operatorname{B}(a,b)$ where
   *  $\operatorname{B}_x(a,b) = \int_0^x t^{a-1} (1-t)^{b-1} dt$ for
   *  $0 \leq x \leq 1$. For $a > 0$, $b > 0$ the value is computed using
   *  a continued fraction expansion. For all other values it is computed
   *  using the relation
   *
   *  .. math::
   *    I_x(a,b) = (1/a) x^a {}_2F_1(a,1-b,a+1,x)/\operatorname{B}(a,b).
   */
  ADDFUNC(gsl_sf_beta_inc, 3);

  /**
   * @file gegenpoly
   *
   * Gegenbauer Functions
   * ====================
   *
   * The Gegenbauer polynomials are defined in Abramowitz & Stegun,
   * Chapter 22, where they are known as Ultraspherical polynomials.
   *
   * .. index::
   *   Gegenbauer function
   *   Gegenbauer polynomial
   */

  /**
   * .. function:: gsl_sf_gegenpoly_1(lambda, x)
   */
  ADDFUNC(gsl_sf_gegenpoly_1, 2);

  /**
   * .. function:: gsl_sf_gegenpoly_2(lambda, x)
   */
  ADDFUNC(gsl_sf_gegenpoly_2, 2);

  /**
   * .. function:: gsl_sf_gegenpoly_3(lambda, x)
   *
   *  These functions evaluate the Gegenbauer polynomials $C^{(\lambda)}_n(x)$
   *  using explicit representations for $n = 1, 2, 3$.
   */
  ADDFUNC(gsl_sf_gegenpoly_3, 2);

  /**
   * .. function:: gsl_sf_gegenpoly_n(n, lambda, x)
   *
   *  This function evaluates the Gegenbauer polynomial $C^{(\lambda)}_n(x)$
   *  for a specific value of $n$, $\lambda$, $x$ subject to $\lambda > -1/2$,
   *  $n \geq 0$.
   */
  ADDFUNC(gsl_sf_gegenpoly_n, 3);

  /**
   * @file hyperg
   *
   * Hypergeometric Functions
   * ========================
   *
   * Hypergeometric functions are described in Abramowitz & Stegun,
   * Chapters 13 and 15.
   *
   * .. index:: Hypergeometric function
   */

  /**
   * .. function:: gsl_sf_hyperg_0F1(c, x)
   *
   *  This routine computes the hypergeometric function ${}_0F_1(c,x)$.
   */
  ADDFUNC(gsl_sf_hyperg_0F1, 2);

  /**
   * .. function:: gsl_sf_hyperg_1F1_int(m, n, x)
   *
   *  This routine computes the confluent hypergeometric function
   *  ${}_1F_1(m,n,x) = M(m,n,x)$ for integer parameters $m$, $n$.
   */
  ADDFUNC(gsl_sf_hyperg_1F1_int, 3);

  /**
   * .. function:: gsl_sf_hyperg_1F1(a, b, x)
   *
   *  This routine computes the confluent hypergeometric function
   *  ${}_1F_1(a,b,x) = M(a,b,x)$ for general parameters $a$, $b$.
   */
  ADDFUNC(gsl_sf_hyperg_1F1, 3);

  /**
   * .. function:: gsl_sf_hyperg_U_int(m, n, x)
   *
   *  This routine computes the confluent hypergeometric function
   *  $U(m,n,x)$ for integer parameters $m$, $n$.
   */
  ADDFUNC(gsl_sf_hyperg_U_int, 3);

  /**
   * .. function:: gsl_sf_hyperg_U(a, b, x)
   *
   *  This routine computes the confluent hypergeometric function $U(a,b,x)$.
   */
  ADDFUNC(gsl_sf_hyperg_U, 3);

  /**
   * .. function:: gsl_sf_hyperg_2F1(a, b, c, x)
   *
   *  This routine computes the Gauss hypergeometric function
   *  ${}_2F_1(a,b,c,x) = F(a,b,c,x)$ for $|x| < 1$.
   *
   *  If the arguments $(a,b,c,x)$ are too close to a singularity then
   *  the function can return an error when the series approximation
   *  converges too slowly. This occurs in the region of
   *  $x=1, c - a - b = m$ for integer $m$.
   */
  ADDFUNC(gsl_sf_hyperg_2F1, 4);

  /**
   * .. function:: gsl_sf_hyperg_2F1_conj(aR, aI, c, x)
   *
   *  This routine computes the Gauss hypergeometric function
   *  ${}_2F_1(a_R + i a_I, a_R - i a_I, c, x)$ with complex parameters
   *  for $|x| < 1$.
   */
  ADDFUNC(gsl_sf_hyperg_2F1_conj, 4);

  /**
   * .. function:: gsl_sf_hyperg_2F1_renorm(a, b, c, x)
   *
   *  This routine computes the renormalized Gauss hypergeometric
   *  function ${}_2F_1(a,b,c,x) / \Gamma(c)$ for $|x| < 1$.
   */
  ADDFUNC(gsl_sf_hyperg_2F1_renorm, 4);

  /**
   * .. function:: gsl_sf_hyperg_2F1_conj_renorm(aR, aI, c, x)
   *
   *  This routine computes the renormalized Gauss hypergeometric
   *  function ${}_2F_1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c)$
   *  for $|x| < 1$.
   */
  ADDFUNC(gsl_sf_hyperg_2F1_conj_renorm, 4);

  /**
   * .. function:: gsl_sf_hyperg_2F0(a, b, x)
   *
   *  This routine computes the hypergeometric function ${}_2F_0(a,b,x)$.
   *  The series representation is a divergent hypergeometric series.
   *  However, for $x < 0$ we have ${}_2F_0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)$
   */
  ADDFUNC(gsl_sf_hyperg_2F0, 3);

  /**
   * @file laguerre
   * .. _laguerre-functions:
   *
   * Laguerre Functions
   * ==================
   *
   * The generalized Laguerre polynomials are defined in terms of confluent
   * hypergeometric functions as $L^a_n(x) = ((a+1)_n / n!) {}_1F_1(-n,a+1,x)$,
   * and are sometimes referred to as the associated Laguerre polynomials.
   * They are related to the plain Laguerre polynomials $L_n(x)$ by
   * $L^0_n(x) = L_n(x)$ and $L^k_n(x) = (-1)^k (d^k/dx^k) L_{n+k}(x)$.
   * For more information see Abramowitz & Stegun, Chapter 22.
   *
   * .. index::
   *   Laguerre function
   *   Laguerre polynomial
   */

  /**
   * .. function:: gsl_sf_laguerre_1(a, x)
   */
  ADDFUNC(gsl_sf_laguerre_1, 2);

  /**
   * .. function:: gsl_sf_laguerre_2(a, x)
   */
  ADDFUNC(gsl_sf_laguerre_2, 2);

  /**
   * .. function:: gsl_sf_laguerre_3(a, x)
   *
   *  These routines evaluate the generalized Laguerre polynomials
   *  $L^a_1(x), L^a_2(x), L^a_3(x)$ using explicit representations.
   */
  ADDFUNC(gsl_sf_laguerre_3, 2);

  /**
   * .. function:: gsl_sf_laguerre_n(n, a, x)
   *
   *  This routine evaluates the generalized Laguerre polynomials
   *  $L^a_n(x)$ for $a > -1, n >= 0$.
   */
  ADDFUNC(gsl_sf_laguerre_n, 3);

  /**
   * @file lambert
   *
   * Lambert W Functions
   * ===================
   *
   * Lambert's $W$ functions, $W(x)$, are defined to be solutions of the
   * equation $W(x) \exp(W(x)) = x$. This function has multiple branches
   * for $x < 0$; however, it has only two real-valued branches. We define
   * $W_0(x)$ to be the principal branch, where $W > -1$ for $x < 0$, and
   * $W_{-1}(x)$ to be the other real branch, where $W < -1$ for $x < 0$.
   *
   * .. index:: Lambert W function
   */

  /**
   * .. function:: gsl_sf_lambert_W0(x)
   *
   *  This routine computes the principal branch of the Lambert $W$ function,
   *  $W_0(x)$.
   */
  ADDFUNC(gsl_sf_lambert_W0, 1);

  /**
   * .. function:: gsl_sf_lambert_Wm1(x)
   *
   *  This routine computes the secondary real-valued branch of the Lambert
   *  $W$ function, $W_{-1}(x)$.
   */
  ADDFUNC(gsl_sf_lambert_Wm1, 1);

  /**
   * @file legendre
   *
   * Legendre Functions and Spherical Harmonics
   * ==========================================
   *
   * The Legendre Functions and Legendre Polynomials are described in
   * Abramowitz & Stegun, Chapter 8.
   *
   * .. index::
   *   Legendre function
   *   Legendre polynomial
   */

  /**
   * Legendre Polynomials
   * --------------------
   */

  /**
   * .. function:: gsl_sf_legendre_P1(x)
   */
  ADDFUNC(gsl_sf_legendre_P1, 1);

  /**
   * .. function:: gsl_sf_legendre_P2(x)
   */
  ADDFUNC(gsl_sf_legendre_P2, 1);

  /**
   * .. function:: gsl_sf_legendre_P3(x)
   *
   *  These functions evaluate the Legendre polynomials $P_l(x)$ using
   *  explicit representations for $l=1, 2, 3$.
   */
  ADDFUNC(gsl_sf_legendre_P3, 1);

  /**
   * .. function:: gsl_sf_legendre_Pl(l, x)
   *
   *  This function evaluates the Legendre polynomial $P_l(x)$ for a
   *  specific value of integer parameter $l$, $x$ subject to
   *  $l \geq 0, |x| \leq 1$.
   */
  ADDFUNC(gsl_sf_legendre_Pl, 2);

  /**
   * .. function:: gsl_sf_legendre_Q0(x)
   *
   *  This routine computes the Legendre function $Q_0(x)$ for
   *  $x > -1, x \ne 1$.
   */
  ADDFUNC(gsl_sf_legendre_Q0, 1);

  /**
   * .. function:: gsl_sf_legendre_Q1(x)
   *
   *  This routine computes the Legendre function $Q_1(x)$ for
   *  $x > -1, x \ne 1$.
   */
  ADDFUNC(gsl_sf_legendre_Q1, 1);

  /**
   * .. function:: gsl_sf_legendre_Ql(l, x)
   *
   *  This routine computes the Legendre function $Q_l(x)$ for
   *  $x > -1, x \ne 1$ and $l \geq 0$.
   */
  ADDFUNC(gsl_sf_legendre_Ql, 2);

  /**
   * Associated Legendre Polynomials and Spherical Harmonics
   * -------------------------------------------------------
   *
   * The following functions compute the associated Legendre Polynomials
   * $P_l^m(x)$. Note that this function grows combinatorially with $l$ and
   * can overflow for $l$ larger than about 150. There is no trouble for
   * small $m$, but overflow occurs when $m$ and $l$ are both large.
   * Rather than allow overflows, these functions refuse to calculate
   * $P_l^m(x)$ and return an error when they can sense that $l$ and $m$
   * are too big.
   *
   * If you want to calculate a spherical harmonic, then do not use these
   * functions. Instead use ``gsl_sf_legendre_sphPlm`` below, which uses a
   * similar recursion, but with the normalized functions.
   */

  /**
   * .. function:: gsl_sf_legendre_Plm(l, m, x)
   *
   *  This routine computes the associated Legendre polynomial
   *  $P_l^m(x)$ for $m \geq 0, l \geq m, |x| \leq 1$.
   */
  ADDFUNC(gsl_sf_legendre_Plm, 3);

  /**
   * .. function:: gsl_sf_legendre_sphPlm(l, m, x)
   *
   *  This routine computes the normalized associated Legendre polynomial
   *  $\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x)$ suitable for use
   *  in spherical harmonics. The parameters must satisfy
   *  $m \geq 0, l \geq m, |x| \leq 1$. Theses routines avoid the overflows
   *  that occur for the standard normalization of $P_l^m(x)$.
   */
  ADDFUNC(gsl_sf_legendre_sphPlm, 3);

  /**
   * Conical Functions
   * -----------------
   *
   * The Conical Functions $P^\mu_{-(1/2)+i\lambda}(x)$ and
   * $Q^\mu_{-(1/2)+i\lambda}$ are described in Abramowitz & Stegun,
   * Section 8.12.
   */

  /**
   * .. function:: gsl_sf_conicalP_half(lambda, x)
   *
   *  This routine computes the irregular Spherical Conical Function
   *  $P^{1/2}_{-1/2 + i \lambda}(x)$ for $x > -1$.
   */
  ADDFUNC(gsl_sf_conicalP_half, 2);

  /**
   * .. function:: gsl_sf_conicalP_mhalf(lambda, x)
   *
   *  This routine computes the regular Spherical Conical Function
   *  $P^{-1/2}_{-1/2 + i \lambda}(x)$ for $x > -1$.
   */
  ADDFUNC(gsl_sf_conicalP_mhalf, 2);

  /**
   * .. function:: gsl_sf_conicalP_0(lambda, x)
   *
   *  This routine computes the conical function
   *  $P^0_{-1/2 + i \lambda}(x)$ for $x > -1$.
   */
  ADDFUNC(gsl_sf_conicalP_0, 2);

  /**
   * .. function:: gsl_sf_conicalP_1(lambda, x)
   *
   *  This routine computes the conical function
   *  $P^1_{-1/2 + i \lambda}(x)$ for $x > -1$.
   */
  ADDFUNC(gsl_sf_conicalP_1, 2);

  /**
   * .. function:: gsl_sf_conicalP_sph_reg(l, lambda, x)
   *
   *  This routine computes the Regular Spherical Conical Function
   *  $P^{-1/2-l}_{-1/2 + i \lambda}(x)$ for $x > -1, l \geq -1$.
   */
  ADDFUNC(gsl_sf_conicalP_sph_reg, 3);

  /**
   * .. function:: gsl_sf_conicalP_cyl_reg(m, lambda, x)
   *
   *  This routine computes the Regular Cylindrical Conical Function
   *  $P^{-m}_{-1/2 + i \lambda}(x)$ for $x > -1, m \geq -1$.
   */
  ADDFUNC(gsl_sf_conicalP_cyl_reg, 3);

  /**
   * Radial Functions for Hyperbolic Space
   * -------------------------------------
   *
   * The following spherical functions are specializations of Legendre
   * functions which give the regular eigenfunctions of the Laplacian
   * on a 3-dimensional hyperbolic space H3d. Of particular interest is
   * the flat limit, $\lambda \to \infty, \eta \to 0, \lambda\eta$ fixed.
   */

  /**
   * .. function:: gsl_sf_legendre_H3d_0(lambda, eta)
   *
   *  This routine computes the zeroth radial eigenfunction of the
   *  Laplacian on the 3-dimensional hyperbolic space,
   *
   *  .. math::
   *    L^{H3d}_0(\lambda,\eta) :=
   *      \sin(\lambda\eta)/(\lambda\sinh(\eta)) \text{ for } \eta \geq 0.
   *
   *  In the flat limit this takes the form
   *  $L^{H3d}_0(\lambda,\eta) = j_0(\lambda\eta)$.
   */
  ADDFUNC(gsl_sf_legendre_H3d_0, 2);

  /**
   * .. function:: gsl_sf_legendre_H3d_1(lambda, eta)
   *
   *  This routine computes the first radial eigenfunction of the
   *  Laplacian on the 3-dimensional hyperbolic space,
   *
   *  .. math::
   *    L^{H3d}_1(\lambda,\eta) := 1/\sqrt{\lambda^2 + 1} \sin(\lambda \eta)/
   *      (\lambda \sinh(\eta)) (\coth(\eta) - \lambda \cot(\lambda\eta))
   *      \text{ for } \eta \geq 0.
   *
   *  In the flat limit this takes the form
   *  $L^{H3d}_1(\lambda,\eta) = j_1(\lambda\eta)$.
   */
  ADDFUNC(gsl_sf_legendre_H3d_1, 2);

  /**
   * .. function:: gsl_sf_legendre_H3d(l, lambda, eta)
   *
   *  This routine computes the $l$-th radial eigenfunction of the
   *  Laplacian on the 3-dimensional hyperbolic space $\eta \geq 0, l \geq 0$.
   *  In the flat limit this takes the form
   *  $L^{H3d}_l(\lambda,\eta) = j_l(\lambda\eta)$.
   */
  ADDFUNC(gsl_sf_legendre_H3d, 3);

  /**
   * @file log
   *
   * Logarithm and Related Functions
   * ===============================
   *
   * Information on the properties of the Logarithm function can be found
   * in Abramowitz & Stegun, Chapter 4.
   *
   * .. index:: logarithm
   */

  /**
   * .. function:: gsl_sf_log(x)
   *
   *  This routine computes the logarithm of $x$, $\log(x)$, for $x > 0$.
   */
  ADDFUNC(gsl_sf_log, 1);

  /**
   * .. function:: gsl_sf_log_abs(x)
   *
   *  This routine computes the logarithm of the magnitude of $x$,
   *  $\log(|x|)$, for $x \ne 0$.
   */
  ADDFUNC(gsl_sf_log_abs, 1);

  /**
   * .. function:: gsl_sf_log_1plusx(x)
   *
   *  This routine computes $\log(1 + x)$ for $x > -1$ using an algorithm
   *  that is accurate for small $x$.
   */
  ADDFUNC(gsl_sf_log_1plusx, 1);

  /**
   * .. function:: gsl_sf_log_1plusx_mx(x)
   *
   *  This routine computes $\log(1 + x) - x$ for $x > -1$ using an
   *  algorithm that is accurate for small $x$.
   */
  ADDFUNC(gsl_sf_log_1plusx_mx, 1);

  /**
   * @file mathieu
   *
   * Mathieu Functions
   * =================
   *
   * The routines described in this section compute the angular and radial
   * Mathieu functions, and their characteristic values. Mathieu functions
   * are the solutions of the following two differential equations:
   *
   * .. math::
   *   d^2y/dv^2 + (a - 2q\cos 2v)y = 0 \\
   *   d^2f/du^2 - (a - 2q\cosh 2u)f = 0
   *
   * .. index:: Mathieu function
   *
   * The angular Mathieu functions $ce_r(x,q), se_r(x,q)$ are the even
   * and odd periodic solutions of the first equation, which is known as
   * Mathieu's equation. These exist only for the discrete sequence of
   * characteristic values $a=a_r(q)$ (even-periodic) and $a=b_r(q)$
   * (odd-periodic).
   *
   * The radial Mathieu functions $Mc^{(j)}_{r}(z,q), Ms^{(j)}_{r}(z,q)$
   * are the solutions of the second equation, which is referred to as
   * Mathieu's modified equation. The radial Mathieu functions of the
   * first, second, third and fourth kind are denoted by the parameter
   * $j$, which takes the value 1, 2, 3 or 4.
   *
   * For more information on the Mathieu functions, see Abramowitz and
   * Stegun, Chapter 20.
   */

  /*
   * Mathieu Function Characteristic Values
   * --------------------------------------
   */

  /**
   * .. function:: gsl_sf_mathieu_a(n, q)
   */
  ADDFUNC(gsl_sf_mathieu_a, 2);

  /**
   * .. function:: gsl_sf_mathieu_b(n, q)
   *
   *  These routines compute the characteristic values $a_n(q), b_n(q)$
   *  of the Mathieu functions $ce_n(q,x)$ and $se_n(q,x)$, respectively.
   */
  ADDFUNC(gsl_sf_mathieu_b, 2);

  /**
   * Angular Mathieu Functions
   * -------------------------
   */

  /**
   * .. function:: gsl_sf_mathieu_ce(n, q, x)
   */
  ADDFUNC(gsl_sf_mathieu_ce, 3);

  /**
   * .. function:: gsl_sf_mathieu_se(n, q, x)
   *
   *  These routines compute the angular Mathieu functions $ce_n(q,x)$ and
   *  $se_n(q,x)$, respectively.
   */
  ADDFUNC(gsl_sf_mathieu_se, 3);

  /**
   * Radial Mathieu Functions
   * ------------------------
   */

  /**
   * .. function:: gsl_sf_mathieu_Mc(j, n, q, x)
   */
  ADDFUNC(gsl_sf_mathieu_Mc, 4);

  /**
   * .. function:: gsl_sf_mathieu_Ms(j, n, q, x)
   *
   *  These routines compute the radial $j$-th kind Mathieu functions
   *  $Mc_n^{(j)}(q,x)$ and $Ms_n^{(j)}(q,x)$ of order $n$.
   *
   *  The allowed values of $j$ are $1$ and $2$. The functions for
   *  $j = 3, 4$ can be computed as $M_n^{(3)} = M_n^{(1)} + iM_n^{(2)}$
   *  and $M_n^{(4)} = M_n^{(1)} - iM_n^{(2)}$, where
   *  $M_n^{(j)} = Mc_n^{(j)}$ or $Ms_n^{(j)}$.
   */
  ADDFUNC(gsl_sf_mathieu_Ms, 4);

  /**
   * @file pow
   *
   * Power Function
   * ==============
   *
   * .. index:: power function
   */

  /**
   * .. function:: gsl_sf_pow_int(x, n)
   *
   *  This routine computes the power $x^n$ for integer $n$. The power is
   *  computed using the minimum number of multiplications. For example,
   *  $x^8$ is computed as $((x^2)^2)^2$, requiring only 3 multiplications.
   *  For reasons of efficiency, these functions do not check for overflow
   *  or underflow conditions.
   *
   *  .. code-block:: none
   *
   *    include gsl.ampl;
   *    # compute and print 3**12
   *    print gsl_sf_pow_int(3, 12);
   */
  ADDFUNC(gsl_sf_pow_int, 2);

  /**
   * @file psi
   *
   * Psi (Digamma) Function
   * ======================
   *
   * The polygamma functions of order $n$ are defined by
   *
   * .. math::
   *   \psi^{(n)}(x) = (d/dx)^n \psi(x) = (d/dx)^{n+1} \log(\Gamma(x))
   *
   * where $\psi(x) = \Gamma'(x)/\Gamma(x)$ is known as the digamma function.
   *
   * .. index::
   *   psi function
   *   polygamma function
   */

  /**
   * Digamma Function
   * ----------------
   */

  /**
   * .. function:: gsl_sf_psi_int(n)
   *
   *  This routine computes the :index:`digamma function` $\psi(n)$ for positive
   *  integer $n$. The digamma function is also called the Psi function.
   */
  ADDFUNC(gsl_sf_psi_int, 1);

  /**
   * .. function:: gsl_sf_psi(x)
   *
   *  This routine computes the digamma function $\psi(x)$ for general
   *  $x, x \ne 0$.
   */
  ADDFUNC(gsl_sf_psi, 1);

  /**
   * .. function:: gsl_sf_psi_1piy(x)
   *
   *  This routine computes the real part of the digamma function on
   *  the line $1+i y, \operatorname{Re}[\psi(1 + i y)]$.
   */
  ADDFUNC(gsl_sf_psi_1piy, 1);

  /**
   * Trigamma Function
   * -----------------
   */

  /**
   * .. function:: gsl_sf_psi_1_int(n)
   *
   *  This routine computes the :index:`Trigamma function` $\psi'(n)$ for
   *  positive integer $n$.
   */
  ADDFUNC(gsl_sf_psi_1_int, 1);

  /**
   * .. function:: gsl_sf_psi_1(x)
   *
   *  This routine computes the Trigamma function $\psi'(x)$ for general $x$.
   */
  ADDFUNC(gsl_sf_psi_1, 1);

  /**
   * Polygamma Function
   * ------------------
   */

  /**
   * .. function:: gsl_sf_psi_n(n, x)
   *
   *  This routine computes the polygamma function $\psi^{(n)}(x)$ for
   *  $n \geq 0, x > 0$.
   */
  ADDFUNC(gsl_sf_psi_n, 2);

  /**
   * @file synchrotron
   *
   * Synchrotron Functions
   * =====================
   */

  /**
   * .. function:: gsl_sf_synchrotron_1(x)
   *
   *  This routine computes the first :index:`synchrotron function`
   *  $x \int_x^\infty K_{5/3}(t) dt$ for $x \geq 0$.
   */
  ADDFUNC(gsl_sf_synchrotron_1, 1);

  /**
   * .. function:: gsl_sf_synchrotron_2(x)
   *
   *  This routine computes the second synchrotron function
   *  $x K_{2/3}(x)$ for $x \geq 0$.
   */
  ADDFUNC(gsl_sf_synchrotron_2, 1);

  /**
   * @file transport
   *
   * Transport Functions
   * ===================
   *
   * The transport functions $J(n,x)$ are defined by the integral
   * representations $J(n,x) := \int_0^x t^n e^t /(e^t - 1)^2 dt$.
   *
   * .. index:: transport function
   */

  /**
   * .. function:: gsl_sf_transport_2(x)
   *
   *  This routine computes the transport function $J(2,x)$.
   */
  ADDFUNC(gsl_sf_transport_2, 1);

  /**
   * .. function:: gsl_sf_transport_3(x)
   *
   *  This routine computes the transport function $J(3,x)$.
   */
  ADDFUNC(gsl_sf_transport_3, 1);

  /**
   * .. function:: gsl_sf_transport_4(x)
   *
   *  This routine computes the transport function $J(4,x)$.
   */
  ADDFUNC(gsl_sf_transport_4, 1);

  /**
   * .. function:: gsl_sf_transport_5(x)
   *
   *  This routine computes the transport function $J(5,x)$.
   */
  ADDFUNC(gsl_sf_transport_5, 1);

  

  /**
  * @file trig
  *
  * Trigonometric Functions
  * =======================
  *
  * A subset of the trigonometric functions defined in GSL is exported here
  *
  * .. index:: trigonometric functions
  */

  /**
 * .. function:: gsl_sf_sinc(x)
 *
 *  This routine computes the :index:`Sinc function`
 *
 *  .. math::
 *    \operatorname{Si}(x) = \int_0^x \sin(t)/t dt.
 */
  ADDFUNC(gsl_sf_sinc, 1);


  /**
   * @file zeta
   *
   * Zeta Functions
   * ==============
   *
   * The :index:`Riemann zeta function` is defined in Abramowitz & Stegun,
   * Section 23.2.
   *
   * .. index:: zeta function
   */

  /**
   * Riemann Zeta Function
   * ---------------------
   *
   * The Riemann zeta function is defined by the infinite sum
   * $\zeta(s) = \sum_{k=1}^\infty k^{-s}$.
   */

  /**
   * .. function:: gsl_sf_zeta_int(n)
   *
   *  This routine computes the Riemann zeta function $\zeta(n)$ for integer
   *  $n, n \ne 1$.
   */
  ADDFUNC(gsl_sf_zeta_int, 1);

  /**
   * .. function:: gsl_sf_zeta(s)
   *
   *  This routine computes the Riemann zeta function $\zeta(s)$ for arbitrary
   *  $s, s \ne 1$.
   */
  ADDFUNC(gsl_sf_zeta, 1);

  /**
   * Riemann Zeta Function Minus One
   * -------------------------------
   *
   * For large positive argument, the Riemann zeta function approaches one.
   * In this region the fractional part is interesting, and therefore we need
   * a function to evaluate it explicitly.
   */

  /**
   * .. function:: gsl_sf_zetam1_int(n)
   *
   *  This routine computes $\zeta(n) - 1$ for integer $n, n \ne 1$.
   */
  ADDFUNC(gsl_sf_zetam1_int, 1);

  /**
   * .. function:: gsl_sf_zetam1(s)
   *
   *  This routine computes $\zeta(s) - 1$ for arbitrary $s, s \ne 1.$.
   */
  ADDFUNC(gsl_sf_zetam1, 1);

  /**
   * Hurwitz Zeta Function
   * ---------------------
   *
   * The :index:`Hurwitz zeta function` is defined by
   * $\zeta(s,q) = \sum_0^\infty (k+q)^{-s}$.
   */

  /**
   * .. function:: gsl_sf_hzeta(s, q)
   *
   *  This routine computes the Hurwitz zeta function $\zeta(s,q)$ for
   *  $s > 1, q > 0$.
   */
  ADDFUNC(gsl_sf_hzeta, 2);

  /**
   * Eta Function
   * ------------
   *
   * The :index:`eta function` is defined by $\eta(s) = (1-2^{1-s}) \zeta(s)$.
   */

  /**
   * .. function:: gsl_sf_eta_int(n)
   *
   *  This routine computes the eta function $\eta(n)$ for integer $n$.
   */
  ADDFUNC(gsl_sf_eta_int, 1);

  /**
   * .. function:: gsl_sf_eta(s)
   *
   *  This routine computes the eta function $\eta(s)$ for arbitrary $s$.
   */
  ADDFUNC(gsl_sf_eta, 1);

  /**
   * @file sf-refs
   *
   * References and Further Reading
   * ------------------------------
   *
   * The library follows the conventions of Abramowitz & Stegun where
   * possible,
   *
   * * Abramowitz & Stegun (eds.), *Handbook of Mathematical Functions*
   *
   * The following papers contain information on the algorithms used to
   * compute the special functions,
   *
   * * Allan J. MacLeod, MISCFUN: A software package to compute uncommon
   *   special functions. *ACM Trans. Math. Soft.*, vol. 22, 1996, 288-301
   * * G.N. Watson, A Treatise on the Theory of Bessel Functions,
   *   2nd Edition (Cambridge University Press, 1944).
   * * G. Nemeth, Mathematical Approximations of Special Functions,
   *   Nova Science Publishers, ISBN 1-56072-052-2
   * * B.C. Carlson, Special Functions of Applied Mathematics (1977)
   * * N. M. Temme, Special Functions: An Introduction to the Classical
   *   Functions of Mathematical Physics (1996), ISBN 978-0471113133.
   * * W.J. Thompson, Atlas for Computing Mathematical Functions, John Wiley
   *   & Sons, New York (1997).
   * * Y.Y. Luke, Algorithms for the Computation of Mathematical Functions,
   *   Academic Press, New York (1977).
   */

  /**
   * @file randist
   *
   * Random Number Distributions
   * ===========================
   *
   * This chapter describes functions for generating random variates and
   * computing their probability distributions. Samples from the distributions
   * described in this chapter can be obtained using any of the random number
   * generators in the library as an underlying source of randomness.
   *
   * .. index::
   *   probability distribution
   *   random number distribution
   *   random number generator
   *
   * In the simplest cases a non-uniform distribution can be obtained
   * analytically from the uniform distribution of a random number generator
   * by applying an appropriate transformation. This method uses one call
   * to the random number generator. More complicated distributions are
   * created by the acceptance-rejection method, which compares the desired
   * distribution against a distribution which is similar and known
   * analytically. This usually requires several samples from the generator.
   *
   * The library also provides cumulative distribution functions and inverse
   * cumulative distribution functions, sometimes referred to as quantile
   * functions. The cumulative distribution functions and their inverses are
   * computed separately for the upper and lower tails of the distribution,
   * allowing full accuracy to be retained for small results.
   *
   * .. index:: cumulative distribution function
   *
   * Note that the discrete random variate functions always return a value
   * of type unsigned int, and on most platforms this has a maximum value of
   * $2^{32}-1 \approx 4.29e9$. They should only be called with a safe range of
   * parameters (where there is a negligible probability of a variate
   * exceeding this limit) to prevent incorrect results due to overflow.
   *
   * .. toctree::
   *    :maxdepth: 2
   *
   *    ran-intro
   *    ran-gaussian
   *    ran-gaussian-tail
   *    ran-exponential
   *    ran-laplace
   *    ran-exppow
   *    ran-cauchy
   *    ran-rayleigh
   *    ran-rayleigh-tail
   *    ran-landau
   *    ran-levy
   *    ran-levy-skew
   *    ran-gamma
   *    ran-flat
   *    ran-lognormal
   *    ran-chisq
   *    ran-fdist
   *    ran-tdist
   *    ran-beta
   *    ran-logistic
   *    ran-pareto
   *    ran-weibull
   *    ran-gumbel1
   *    ran-gumbel2
   *    ran-poisson
   *    ran-bernoulli
   *    ran-binomial
   *    ran-negative-binomial
   *    ran-pascal
   *    ran-geometric
   *    ran-hypergeometric
   *    ran-logarithmic
   *    ran-refs
   */

  /**
   * @file ran-intro
   *
   * Introduction
   * ============
   *
   * Continuous random number distributions are defined by a probability
   * density function, $p(x)$, such that the probability of $x$ occurring
   * in the infinitesimal range $x$ to $x+dx$ is $p dx$.
   *
   * The cumulative distribution function for the lower tail $P(x)$ is defined
   * by the integral,
   *
   * .. math::
   *   P(x) = \int_{-\infty}^{x} dx' p(x')
   *
   * and gives the probability of a variate taking a value less than $x$.
   *
   * The cumulative distribution function for the upper tail $Q(x)$ is defined
   * by the integral,
   *
   * .. math::
   *   Q(x) = \int_{x}^{+\infty} dx' p(x')
   *
   * and gives the probability of a variate taking a value greater than $x$.
   *
   * The upper and lower cumulative distribution functions are related by
   * $P(x) + Q(x) = 1$ and satisfy $0 \leq P(x) \leq 1, 0 \leq Q(x) \leq 1$.
   *
   * The inverse cumulative distributions, $x=P^{-1}(P)$ and $x=Q^{-1}(Q)$
   * give the values of $x$ which correspond to a specific value of $P$ or $Q$.
   * They can be used to find confidence limits from probability values.
   *
   * For discrete distributions the probability of sampling the integer
   * value $k$ is given by $p(k)$, where $\sum_k p(k) = 1$. The cumulative
   * distribution for the lower tail $P(k)$ of a discrete distribution is
   * defined as,
   *
   * .. math::
   *   P(k) = \sum_{i \leq k} p(i)
   *
   * where the sum is over the allowed range of the distribution less than
   * or equal to $k$.
   *
   * The cumulative distribution for the upper tail of a discrete distribution
   * $Q(k)$ is defined as
   *
   * .. math::
   *   Q(k) = \sum_{i > k} p(i)
   *
   * giving the sum of probabilities for all values greater than $k$.
   * These two definitions satisfy the identity $P(k)+Q(k)=1$.
   *
   * If the range of the distribution is $1$ to $n$ inclusive then
   * $P(n)=1, Q(n)=0$ while $P(1) = p(1), Q(1)=1-p(1)$.
   */

  /* Initialize the random number generator. */
#ifdef addrandinit
  if (ae->ASLdate >= 20120830)
  	addrandinit(rng_init, ae);
  else
#endif
  rng = gsl_rng_alloc(gsl_rng_env_setup());
  at_reset(free_rng, &rng);

  /**
   * @file ran-gaussian
   *
   * The Gaussian Distribution
   * =========================
   *
   * .. index:: Gaussian distribution
   */

  /**
   * .. function:: gsl_ran_gaussian(sigma)
   *
   *  This function returns a :index:`Gaussian random variate`, with mean
   *  zero and standard deviation ``sigma``. The probability distribution
   *  for Gaussian random variates is,
   *
   *  .. math::
   *    p(x) dx = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2) dx
   *
   *  for $x$ in the range $-\infty$ to $+\infty$. Use the transformation
   *  $z = \mu + x$ on the numbers returned by ``gsl_ran_gaussian`` to obtain
   *  a Gaussian distribution with mean $\mu$. This function uses the
   *  Box-Muller algorithm which requires two calls to the random number
   *  generator.
   */
  ADDFUNC_RANDOM(gsl_ran_gaussian, 1);

  /**
   * .. function:: gsl_ran_gaussian_pdf(x, sigma)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Gaussian distribution with standard deviation ``sigma``, using the formula
   *  given above.
   */
  ADDFUNC(gsl_ran_gaussian_pdf, 2);

  /**
   * .. function:: gsl_ran_gaussian_ziggurat(sigma)
   */
  ADDFUNC_RANDOM(gsl_ran_gaussian_ziggurat, 1);

  /**
   * .. function:: gsl_ran_gaussian_ratio_method(sigma)
   *
   *  These functions compute a Gaussian random variate using the alternative
   *  Marsaglia-Tsang ziggurat and Kinderman-Monahan-Leva ratio methods.
   *  The Ziggurat algorithm is the fastest available algorithm in most cases.
   */
  ADDFUNC_RANDOM(gsl_ran_gaussian_ratio_method, 1);

  /**
   * .. function:: gsl_ran_ugaussian()
   */
  ADDFUNC_RANDOM(gsl_ran_ugaussian, 0);

  /**
   * .. function:: gsl_ran_ugaussian_pdf(x)
   */
  ADDFUNC(gsl_ran_ugaussian_pdf, 1);

  /**
   * .. function:: gsl_ran_ugaussian_ratio_method()
   *
   *  These functions compute results for the unit Gaussian distribution.
   *  They are equivalent to the functions above with a standard deviation
   *  of one, ``sigma`` = 1.
   */
  ADDFUNC_RANDOM(gsl_ran_ugaussian_ratio_method, 0);

  /**
   * .. function:: gsl_cdf_gaussian_P(x, sigma)
   */
  ADDFUNC(gsl_cdf_gaussian_P, 2);

  /**
   * .. function:: gsl_cdf_gaussian_Q(x, sigma)
   */
  ADDFUNC(gsl_cdf_gaussian_Q, 2);

  /**
   * .. function:: gsl_cdf_gaussian_Pinv(P, sigma)
   */
  ADDFUNC(gsl_cdf_gaussian_Pinv, 2);

  /**
   * .. function:: gsl_cdf_gaussian_Qinv(Q, sigma)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the Gaussian distribution with
   *  standard deviation ``sigma``.
   */
  ADDFUNC(gsl_cdf_gaussian_Qinv, 2);

  /**
   * .. function:: gsl_cdf_ugaussian_P(x)
   */
  ADDFUNC(gsl_cdf_ugaussian_P, 1);

  /**
   * .. function:: gsl_cdf_ugaussian_Q(x)
   */
  ADDFUNC(gsl_cdf_ugaussian_Q, 1);

  /**
   * .. function:: gsl_cdf_ugaussian_Pinv(P)
   */
  ADDFUNC(gsl_cdf_ugaussian_Pinv, 1);

  /**
   * .. function:: gsl_cdf_ugaussian_Qinv(Q)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the unit Gaussian distribution.
   */
  ADDFUNC(gsl_cdf_ugaussian_Qinv, 1);

  /**
   * @file ran-gaussian-tail
   *
   * The Gaussian Tail Distribution
   * ==============================
   */

  /**
   * .. function:: gsl_ran_gaussian_tail(a, sigma)
   *
   *  This function provides random variates from the upper tail of a
   *  Gaussian distribution with standard deviation ``sigma``. The values
   *  returned are larger than the lower limit ``a``, which must be positive.
   *  The method is based on Marsaglia's famous rectangle-wedge-tail
   *  algorithm (Ann. Math. Stat. 32, 894-899 (1961)), with this aspect
   *  explained in Knuth, v2, 3rd ed, p139,586 (exercise 11).
   *
   *  The probability distribution for Gaussian tail random variates is,
   *
   *  .. math::
   *    p(x) dx = {1 \over N(a;\sigma) \sqrt{2 \pi \sigma^2}}
   *      \exp (- x^2/(2 \sigma^2)) dx
   *
   *  for $x > a$ where $N(a;\sigma)$ is the normalization constant,
   *
   *  .. math::
   *    N(a;\sigma) = (1/2) \operatorname{erfc}(a / \sqrt{2 \sigma^2}).
   */
  ADDFUNC_RANDOM(gsl_ran_gaussian_tail, 2);

  /**
   * .. function:: gsl_ran_gaussian_tail_pdf(x, a, sigma)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Gaussian tail distribution with standard deviation ``sigma`` and lower
   *  limit ``a``, using the formula given above.
   */
  ADDFUNC(gsl_ran_gaussian_tail_pdf, 3);

  /**
   * .. function:: gsl_ran_ugaussian_tail(a)
   */
  ADDFUNC_RANDOM(gsl_ran_ugaussian_tail, 1);

  /**
   * .. function:: gsl_ran_ugaussian_tail_pdf(x, a)
   *
   *  These functions compute results for the tail of a unit Gaussian
   *  distribution. They are equivalent to the functions above with a
   *  standard deviation of one, ``sigma`` = 1.
   */
  ADDFUNC(gsl_ran_ugaussian_tail_pdf, 2);

  /* The bivariate Gaussian distribution is not wrapped because it returns
     more than one value. */

  /**
   * @file ran-exponential
   *
   * The Exponential Distribution
   * ============================
   */

  /**
   * .. function:: gsl_ran_exponential(mu)
   *
   *  This function returns a random variate from the :index:`exponential
   *  distribution` with mean ``mu``. The distribution is,
   *
   *  .. math::
   *    p(x) dx = {1 \over \mu} \exp(-x/\mu) dx
   *
   *  for $x \geq 0$.
   */
  ADDFUNC_RANDOM(gsl_ran_exponential, 1);

  /**
   * .. function:: gsl_ran_exponential_pdf(x, mu)
   *
   *  This function computes the probability density $p(x)$ at $x$ for an
   *  exponential distribution with mean ``mu``, using the formula given above.
   */
  ADDFUNC(gsl_ran_exponential_pdf, 2);

  /**
   * .. function:: gsl_cdf_exponential_P(x, mu)
   */
  ADDFUNC(gsl_cdf_exponential_P, 2);

  /**
   * .. function:: gsl_cdf_exponential_Q(x, mu)
   */
  ADDFUNC(gsl_cdf_exponential_Q, 2);

  /**
   * .. function:: gsl_cdf_exponential_Pinv(P, mu)
   */
  ADDFUNC(gsl_cdf_exponential_Pinv, 2);

  /**
   * .. function:: gsl_cdf_exponential_Qinv(Q, mu)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the exponential distribution
   *  with mean ``mu``.
   */
  ADDFUNC(gsl_cdf_exponential_Qinv, 2);

  /**
   * @file ran-laplace
   *
   * The Laplace Distribution
   * ========================
   */

  /**
   * .. function:: gsl_ran_laplace(a)
   *
   *  This function returns a random variate from the :index:`Laplace
   *  distribution` with width ``a``. The distribution is,
   *
   *  .. math::
   *    p(x) dx = {1 \over 2 a}  \exp(-|x/a|) dx
   *
   *  for $-\infty < x < \infty$.
   */
  ADDFUNC_RANDOM(gsl_ran_laplace, 1);

  /**
   * .. function:: gsl_ran_laplace_pdf(x, a)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Laplace distribution with width ``a``, using the formula given above.
   */
  ADDFUNC(gsl_ran_laplace_pdf, 2);

  /**
   * .. function:: gsl_ran_laplace_P(x, a)
   */
  ADDFUNC(gsl_cdf_laplace_P, 2);

  /**
   * .. function:: gsl_ran_laplace_Q(x, a)
   */
  ADDFUNC(gsl_cdf_laplace_Q, 2);

  /**
   * .. function:: gsl_ran_laplace_Pinv(P, a)
   */
  ADDFUNC(gsl_cdf_laplace_Pinv, 2);

  /**
   * .. function:: gsl_ran_laplace_Qinv(Q, a)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the Laplace distribution
   *  with width ``a``.
   */
  ADDFUNC(gsl_cdf_laplace_Qinv, 2);

  /**
   * @file ran-exppow
   *
   * The Exponential Power Distribution
   * ==================================
   */

  /**
   * .. function:: gsl_ran_exppow(a, b)
   *
   *  This function returns a random variate from :index:`exponential power
   *  distribution` with scale parameter ``a`` and exponent ``b``.
   *  The distribution is,
   *
   *  .. math::
   *    p(x) dx = {1 \over 2 a \Gamma(1+1/b)} \exp(-|x/a|^b) dx
   *
   *  for $x \geq 0$. For $b = 1$ this reduces to the Laplace distribution.
   *  For $b = 2$ it has the same form as a Gaussian distribution, but with
   *  $a = \sqrt{2} \sigma$.
   */
  ADDFUNC_RANDOM(gsl_ran_exppow, 2);

  /**
   * .. function:: gsl_ran_exppow_pdf(x, a, b)
   *
   *  This function computes the probability density $p(x)$ at $x$ for an
   *  exponential power distribution with scale parameter ``a`` and exponent
   *  ``b``, using the formula given above.
   */
  ADDFUNC(gsl_ran_exppow_pdf, 3);

  /**
   * .. function:: gsl_ran_exppow_P(x, a, b)
   */
  ADDFUNC(gsl_cdf_exppow_P, 3);

  /**
   * .. function:: gsl_ran_exppow_Q(x, a, b)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ for the exponential power distribution with parameters
   *  ``a`` and ``b``.
   */
  ADDFUNC(gsl_cdf_exppow_Q, 3);

  /**
   * @file ran-cauchy
   *
   * The Cauchy Distribution
   * ========================
   */

  /**
   * .. function:: gsl_ran_cauchy(a)
   *
   *  This function returns a random variate from the :index:`Cauchy
   *  distribution` with scale parameter ``a``. The probability
   *  distribution for Cauchy random variates is,
   *
   *  .. math::
   *    p(x) dx = {1 \over a\pi (1 + (x/a)^2) } dx
   *
   *  for $x$ in the range $-\infty$ to $+\infty$. The Cauchy distribution
   *  is also known as the Lorentz distribution.
   */
  ADDFUNC_RANDOM(gsl_ran_cauchy, 1);

  /**
   * .. function:: gsl_ran_cauchy_pdf(x, a)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Cauchy distribution with scale parameter ``a``, using the formula
   *  given above.
   */
  ADDFUNC(gsl_ran_cauchy_pdf, 2);

  /**
   * .. function:: gsl_ran_cauchy_P(x, a)
   */
  ADDFUNC(gsl_cdf_cauchy_P, 2);

  /**
   * .. function:: gsl_ran_cauchy_Q(x, a)
   */
  ADDFUNC(gsl_cdf_cauchy_Q, 2);

  /**
   * .. function:: gsl_ran_cauchy_Pinv(P, a)
   */
  ADDFUNC(gsl_cdf_cauchy_Pinv, 2);

  /**
   * .. function:: gsl_ran_cauchy_Qinv(Q, a)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the Cauchy distribution with
   *  scale parameter ``a``.
   */
  ADDFUNC(gsl_cdf_cauchy_Qinv, 2);

  /**
   * @file ran-rayleigh
   *
   * The Rayleigh Distribution
   * =========================
   */

  /**
   * .. function:: gsl_ran_rayleigh(sigma)
   *
   *  This function returns a random variate from the :index:`Rayleigh
   *  distribution` with scale parameter ``sigma``.
   *  The distribution is,
   *
   *  .. math::
   *    p(x) dx = {x \over \sigma^2} \exp(- x^2/(2 \sigma^2)) dx
   *
   *  for $x > 0$.
   */
  ADDFUNC_RANDOM(gsl_ran_rayleigh, 1);

  /**
   * .. function:: gsl_ran_rayleigh_pdf(x, sigma)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Rayleigh distribution with scale parameter ``sigma``, using the formula
   *  given above.
   */
  ADDFUNC(gsl_ran_rayleigh_pdf, 2);

  /**
   * .. function:: gsl_ran_rayleigh_P(x, sigma)
   */
  ADDFUNC(gsl_cdf_rayleigh_P, 2);

  /**
   * .. function:: gsl_ran_rayleigh_Q(x, sigma)
   */
  ADDFUNC(gsl_cdf_rayleigh_Q, 2);

  /**
   * .. function:: gsl_ran_rayleigh_Pinv(P, sigma)
   */
  ADDFUNC(gsl_cdf_rayleigh_Pinv, 2);

  /**
   * .. function:: gsl_ran_rayleigh_Qinv(Q, sigma)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the Rayleigh distribution with
   *  scale parameter ``sigma``.
   */
  ADDFUNC(gsl_cdf_rayleigh_Qinv, 2);

  /**
   * @file ran-rayleigh-tail
   *
   * The Rayleigh Tail Distribution
   * ==============================
   */

  /**
   * .. function:: gsl_ran_rayleigh_tail(a, sigma)
   *
   *  This function returns a random variate from the tail of the Rayleigh
   *  distribution with scale parameter ``sigma`` and a lower limit of
   *  ``a``. The distribution is,
   *
   *  .. math::
   *    p(x) dx = {x \over \sigma^2} \exp ((a^2 - x^2) /(2 \sigma^2)) dx
   *
   *  for $x > a$.
   */
  ADDFUNC_RANDOM(gsl_ran_rayleigh_tail, 2);

  /**
   * .. function:: gsl_ran_rayleigh_tail_pdf(x, a, sigma)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Rayleigh tail distribution with scale parameter ``sigma`` and lower
   *  limit ``a``, using the formula given above.
   */
  ADDFUNC(gsl_ran_rayleigh_tail_pdf, 3);

  /**
   * @file ran-landau
   *
   * The Landau Distribution
   * =======================
   */

  /**
   * .. function:: gsl_ran_landau()
   *
   *  This function returns a random variate from the :index:`Landau
   *  distribution`. The probability distribution for Landau random
   *  variates is defined analytically by the complex integral,
   *
   *  .. math::
   *    p(x) = (1/(2 \pi i))
   *      \int_{c-i\infty}^{c+i\infty} \exp(s \log(s) + x s) ds
   *
   *  For numerical purposes it is more convenient to use the following
   *  equivalent form of the integral,
   *
   *  .. math::
   *    p(x) = (1/\pi) \int_0^\infty \exp(-t \log(t) - x t) \sin(\pi t) dt.
   */
  ADDFUNC_RANDOM(gsl_ran_landau, 0);

  /**
   * .. function:: gsl_ran_landau_pdf(x)
   *
   *  This function computes the probability density $p(x)$ at $x$ for the
   *  Landau distribution using an approximation to the formula given above.
   */
  ADDFUNC(gsl_ran_landau_pdf, 1);

  /**
   * @file ran-levy
   *
   * The Levy alpha-Stable Distribution
   * ==================================
   */

  /**
   * .. function:: gsl_ran_levy(c, alpha)
   *
   *  This function returns a random variate from the :index:`Levy symmetric
   *  stable distribution` with scale ``c`` and exponent ``alpha``.
   *  The symmetric stable probability distribution is defined by a
   *  Fourier transform,
   *
   *  .. math::
   *    p(x) = {1 \over 2 \pi}
   *      \int_{-\infty}^{+\infty} \exp(-it x - |c t|^\alpha) dt
   *
   *  There is no explicit solution for the form of $p(x)$ and the library
   *  does not define a corresponding pdf function. For $\alpha = 1$ the
   *  distribution reduces to the Cauchy distribution. For $\alpha = 2$ it
   *  is a Gaussian distribution with $\sigma = \sqrt{2} c$. For $\alpha < 1$
   *  the tails of the distribution become extremely wide.
   *
   *  The algorithm only works for $0 < \alpha \leq 2$.
   */
  ADDFUNC_RANDOM(gsl_ran_levy, 2);

  /**
   * @file ran-levy-skew
   *
   * The Levy skew alpha-Stable Distribution
   * =======================================
   */

  /**
   * .. function:: gsl_ran_levy_skew(c, alpha, beta)
   *
   *  This function returns a random variate from the :index:`Levy skew stable
   *  distribution` with scale ``c``, exponent ``alpha`` and skewness
   *  parameter ``beta``. The skewness parameter must lie in the range
   *  [-1,1]. The Levy skew stable probability distribution is defined
   *  by a Fourier transform,
   *
   *  .. math::
   *    p(x) = {1 \over 2 \pi} \int_{-\infty}^{+\infty}
   *      \exp(-it x - |c t|^\alpha (1-i \beta \operatorname{sign}(t)
   *      \tan(\pi \alpha/2))) dt
   *
   *  When $\alpha = 1$ the term $\tan(\pi \alpha/2)$ is replaced by
   *  $-(2/\pi)\log|t|$. There is no explicit solution for the form of
   *  $p(x)$ and the library does not define a corresponding pdf function.
   *  For $\alpha = 2$ the distribution reduces to a Gaussian distribution
   *  with $\sigma = \sqrt{2} c$ and the skewness parameter has no effect.
   *  For $\alpha < 1$ the tails of the distribution become extremely wide.
   *  The symmetric distribution corresponds to $\beta = 0$.
   *
   *  The algorithm only works for $0 < \alpha \leq 2$.
   *
   * The Levy alpha-stable distributions have the property that if
   * $N$ alpha-stable variates are drawn from the distribution
   * $p(c, \alpha, \beta)$ then the sum $Y = X_1 + X_2 + \dots + X_N$
   * will also be distributed as an alpha-stable variate,
   * $p(N^{1/\alpha} c, \alpha, \beta)$.
   */
  ADDFUNC_RANDOM(gsl_ran_levy_skew, 3);

  /**
   * @file ran-gamma
   *
   * The Gamma Distribution
   * ======================
   */

  /**
   * .. function:: gsl_ran_gamma(a, b)
   *
   *  This function returns a random variate from the :index:`gamma
   *  distribution`. The distribution function is,
   *
   *  .. math::
   *    p(x) dx = {1 \over \Gamma(a) b^a} x^{a-1} e^{-x/b} dx
   *
   *  for $x > 0$.
   *
   *  The gamma distribution with an integer parameter ``a`` is known as
   *  the Erlang distribution.
   *
   *  The variates are computed using the Marsaglia-Tsang fast gamma method.
   */
  ADDFUNC_RANDOM(gsl_ran_gamma, 2);

  /**
   * .. function:: gsl_ran_gamma_knuth(a, b)
   *
   *  This function returns a gamma variate using the algorithms from Knuth
   *  (vol 2).
   */
  ADDFUNC_RANDOM(gsl_ran_gamma_knuth, 2);

  /**
   * .. function:: gsl_ran_gamma_pdf(x, a, b)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  gamma distribution with parameters ``a`` and ``b``, using the formula
   *  given above.
   */
  ADDFUNC(gsl_ran_gamma_pdf, 3);

  /**
   * .. function:: gsl_cdf_gamma_P(x, a, b)
   */
  ADDFUNC(gsl_cdf_gamma_P, 3);

  /**
   * .. function:: gsl_cdf_gamma_Q(x, a, b)
   */
  ADDFUNC(gsl_cdf_gamma_Q, 3);

  /**
   * .. function:: gsl_cdf_gamma_Pinv(P, a, b)
   */
  ADDFUNC(gsl_cdf_gamma_Pinv, 3);

  /**
   * .. function:: gsl_cdf_gamma_Qinv(Q, a, b)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the gamma distribution with
   *  parameters ``a`` and ``b``.
   */
  ADDFUNC(gsl_cdf_gamma_Qinv, 3);

  /**
   * @file ran-flat
   *
   * The Flat (Uniform) Distribution
   * ===============================
   */

  /**
   * .. function:: gsl_ran_flat(a, b)
   *
   *  This function returns a random variate from the flat (uniform)
   *  distribution from ``a`` to ``b``. The distribution is,
   *
   *  .. math::
   *    p(x) dx = {1 \over (b-a)} dx
   *
   *  if $a \leq x < b$ and $0$ otherwise.
   *
   *  .. index::
   *    flat distribution
   *    uniform distribution
   */
  ADDFUNC_RANDOM(gsl_ran_flat, 2);

  /**
   * .. function:: gsl_ran_flat_pdf(x, a, b)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  uniform distribution from ``a`` to ``b``, using the formula given above.
   */
  ADDFUNC(gsl_ran_flat_pdf, 3);

  /**
   * .. function:: gsl_cdf_flat_P(x, a, b)
   */
  ADDFUNC(gsl_cdf_flat_P, 3);

  /**
   * .. function:: gsl_cdf_flat_Q(x, a, b)
   */
  ADDFUNC(gsl_cdf_flat_Q, 3);

  /**
   * .. function:: gsl_cdf_flat_Pinv(P, a, b)
   */
  ADDFUNC(gsl_cdf_flat_Pinv, 3);

  /**
   * .. function:: gsl_cdf_flat_Qinv(Q, a, b)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for a uniform distribution from
   *  ``a`` to ``b``.
   */
  ADDFUNC(gsl_cdf_flat_Qinv, 3);

  /**
   * @file ran-lognormal
   *
   * The Lognormal Distribution
   * ==========================
   */

  /**
   * .. function:: gsl_ran_lognormal(zeta, sigma)
   *
   *  This function returns a random variate from the :index:`lognormal
   *  distribution`. The distribution function is,
   *
   *  .. math::
   *    p(x) dx = {1 \over x \sqrt{2 \pi \sigma^2} }
   *      \exp(-(\ln(x) - \zeta)^2/2 \sigma^2) dx
   *
   *  for $x > 0$.
   */
  ADDFUNC_RANDOM(gsl_ran_lognormal, 2);

  /**
   * .. function:: gsl_ran_lognormal_pdf(x, zeta, sigma)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  lognormal distribution with parameters ``zeta`` and ``sigma``, using
   *  the formula given above.
   */
  ADDFUNC(gsl_ran_lognormal_pdf, 3);

  /**
   * .. function:: gsl_cdf_lognormal_P(x, zeta, sigma)
   */
  ADDFUNC(gsl_cdf_lognormal_P, 3);

  /**
   * .. function:: gsl_cdf_lognormal_Q(x, zeta, sigma)
   */
  ADDFUNC(gsl_cdf_lognormal_Q, 3);

  /**
   * .. function:: gsl_cdf_lognormal_Pinv(P, zeta, sigma)
   */
  ADDFUNC(gsl_cdf_lognormal_Pinv, 3);

  /**
   * .. function:: gsl_cdf_lognormal_Qinv(Q, zeta, sigma)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the lognormal distribution
   *  with parameters ``zeta`` and ``sigma``.
   */
  ADDFUNC(gsl_cdf_lognormal_Qinv, 3);

  /**
   * @file ran-chisq
   *
   * The Chi-squared Distribution
   * ============================
   *
   * The :index:`chi-squared distribution` arises in statistics. If $Y_i$
   * are $n$ independent Gaussian random variates with unit variance then the
   * sum-of-squares,
   *
   * .. math::
   *   X_i = \sum_i Y_i^2
   *
   * has a chi-squared distribution with $n$ degrees of freedom.
   */

  /**
   * .. function:: gsl_ran_chisq(nu)
   *
   *  This function returns a random variate from chi-squared distribution
   *  with ``nu`` degrees of freedom. The distribution function is,
   *
   *  .. math::
   *    p(x) dx = {1 \over 2 \Gamma(\nu/2) } (x/2)^{\nu/2 - 1} \exp(-x/2) dx
   *
   *  for $x \geq 0$.
   */
  ADDFUNC_RANDOM(gsl_ran_chisq, 1);

  /**
   * .. function:: gsl_ran_chisq_pdf(x, nu)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  chi-squared distribution with ``nu`` degrees of freedom, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_chisq_pdf, 2);

  /**
   * .. function:: gsl_ran_chisq_P(x, nu)
   */
  ADDFUNC(gsl_cdf_chisq_P, 2);

  /**
   * .. function:: gsl_ran_chisq_Q(x, nu)
   */
  ADDFUNC(gsl_cdf_chisq_Q, 2);

  /**
   * .. function:: gsl_ran_chisq_Pinv(P, nu)
   */
  ADDFUNC(gsl_cdf_chisq_Pinv, 2);

  /**
   * .. function:: gsl_ran_chisq_Qinv(Q, nu)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the chi-squared distribution
   *  with ``nu`` degrees of freedom.
   */
  ADDFUNC(gsl_cdf_chisq_Qinv, 2);

  /**
   * @file ran-fdist
   *
   * The F-distribution
   * ==================
   *
   * The :index:`F-distribution` arises in statistics. If $Y_1$ and $Y_2$
   * are chi-squared deviates with $\nu_1$ and $\nu_2$ degrees of freedom
   * then the ratio,
   *
   * .. math::
   *   X = { (Y_1 / \nu_1) \over (Y_2 / \nu_2) }
   *
   * has an F-distribution $F(x;\nu_1,\nu_2)$.
   */

  /**
   * .. function:: gsl_ran_fdist(nu1, nu2)
   *
   *  This function returns a random variate from the F-distribution with
   *  degrees of freedom ``nu1`` and ``nu2``. The distribution function is,
   *
   *  .. math::
   *    p(x) dx =
   *         { \Gamma((\nu_1 + \nu_2)/2)
   *              \over \Gamma(\nu_1/2) \Gamma(\nu_2/2) }
   *         \nu_1^{\nu_1/2} \nu_2^{\nu_2/2}
   *         x^{\nu_1/2 - 1} (\nu_2 + \nu_1 x)^{-\nu_1/2 -\nu_2/2}
   *
   * for $x \geq 0$.
   */
  ADDFUNC_RANDOM(gsl_ran_fdist, 2);

  /**
   * .. function:: gsl_ran_fdist_pdf(x, nu1, nu2)
   *
   *  This function computes the probability density $p(x)$ at $x$ for an
   *  F-distribution with ``nu1`` and ``nu2`` degrees of freedom, using
   *  the formula given above.
   */
  ADDFUNC(gsl_ran_fdist_pdf, 3);

  /**
   * .. function:: gsl_cdf_fdist_P(x, nu1, nu2)
   */
  ADDFUNC(gsl_cdf_fdist_P, 3);

  /**
   * .. function:: gsl_cdf_fdist_Q(x, nu1, nu2)
   */
  ADDFUNC(gsl_cdf_fdist_Q, 3);

  /**
   * .. function:: gsl_cdf_fdist_Pinv(P, nu1, nu2)
   */
  ADDFUNC(gsl_cdf_fdist_Pinv, 3);

  /**
   * .. function:: gsl_cdf_fdist_Qinv(Q, nu1, nu2)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the F-distribution with
   *  ``nu1`` and ``nu2`` degrees of freedom.
   */
  ADDFUNC(gsl_cdf_fdist_Qinv, 3);

  /**
   * @file ran-tdist
   *
   * The t-distribution
   * ============================
   *
   * The :index:`t-distribution` arises in statistics. If $Y_1$ has a normal
   * distribution and $Y_2$ has a chi-squared distribution with $\nu$
   * degrees of freedom then the ratio,
   *
   * .. math::
   *   X = { Y_1 \over \sqrt{Y_2 / \nu} }
   *
   * has a t-distribution $t(x;\nu)$ with $\nu$ degrees of freedom.
   */

  /**
   * .. function:: gsl_ran_tdist(nu)
   *
   *  This function returns a random variate from the t-distribution.
   *  The distribution function is,
   *
   *  .. math::
   *    p(x) dx = {\Gamma((\nu + 1)/2) \over \sqrt{\pi \nu} \Gamma(\nu/2)}
   *         (1 + x^2/\nu)^{-(\nu + 1)/2} dx
   *
   *  for $-\infty < x < +\infty$.
   */
  ADDFUNC_RANDOM(gsl_ran_tdist, 1);

  /**
   * .. function:: gsl_ran_tdist_pdf(x, nu)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  t-distribution with ``nu`` degrees of freedom, using the formula
   *  given above.
   */
  ADDFUNC(gsl_ran_tdist_pdf, 2);

  /**
   * .. function:: gsl_cdf_tdist_P(x, nu)
   */
  ADDFUNC(gsl_cdf_tdist_P, 2);

  /**
   * .. function:: gsl_cdf_tdist_Q(x, nu)
   */
  ADDFUNC(gsl_cdf_tdist_Q, 2);

  /**
   * .. function:: gsl_cdf_tdist_Pinv(P, nu)
   */
  ADDFUNC(gsl_cdf_tdist_Pinv, 2);

  /**
   * .. function:: gsl_cdf_tdist_Qinv(Q, nu)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the t-distribution with ``nu``
   *  degrees of freedom.
   */
  ADDFUNC(gsl_cdf_tdist_Qinv, 2);

  /**
   * @file ran-beta
   *
   * The Beta Distribution
   * =====================
   */

  /**
   * .. function:: gsl_ran_beta(a, b)
   *
   *  This function returns a random variate from the :index:`beta
   *  distribution`. The distribution function is,
   *
   *  .. math::
   *    p(x) dx = {\Gamma(a+b) \over \Gamma(a) \Gamma(b)}
   *      x^{a-1} (1-x)^{b-1} dx
   *
   *  for $0 \leq x \leq 1$.
   */
  ADDFUNC_RANDOM(gsl_ran_beta, 2);

  /**
   * .. function:: gsl_ran_beta_pdf(x, a, b)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  beta distribution with parameters ``a`` and ``b``, using the formula
   *  given above.
   */
  ADDFUNC(gsl_ran_beta_pdf, 3);

  /**
   * .. function:: gsl_cdf_beta_P(x, a, b)
   */
  ADDFUNC(gsl_cdf_beta_P, 3);

  /**
   * .. function:: gsl_cdf_beta_Q(x, a, b)
   */
  ADDFUNC(gsl_cdf_beta_Q, 3);

  /**
   * .. function:: gsl_cdf_beta_Pinv(P, a, b)
   */
  ADDFUNC(gsl_cdf_beta_Pinv, 3);

  /**
   * .. function:: gsl_cdf_beta_Qinv(Q, a, b)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the beta distribution with
   *  parameters ``a`` and ``b``.
   */
  ADDFUNC(gsl_cdf_beta_Qinv, 3);

  /**
   * @file ran-logistic
   *
   * The Logistic Distribution
   * =========================
   */

  /**
   * .. function:: gsl_ran_logistic(a)
   *
   *  This function returns a random variate from the :index:`logistic
   *  distribution`. The distribution is,
   *
   *  .. math::
   *    p(x) dx = { \exp(-x/a) \over a (1 + \exp(-x/a))^2 } dx
   *
   *  for $-\infty < x < \infty$.
   */
  ADDFUNC_RANDOM(gsl_ran_logistic, 1);

  /**
   * .. function:: gsl_ran_logistic_pdf(x, a)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  logistic distribution with scale parameter a, using the formula given
   *  above.
   */
  ADDFUNC(gsl_ran_logistic_pdf, 2);

  /**
   * .. function:: gsl_ran_logistic_P(x, a)
   */
  ADDFUNC(gsl_cdf_logistic_P, 2);

  /**
   * .. function:: gsl_ran_logistic_Q(x, a)
   */
  ADDFUNC(gsl_cdf_logistic_Q, 2);

  /**
   * .. function:: gsl_ran_logistic_Pinv(P, a)
   */
  ADDFUNC(gsl_cdf_logistic_Pinv, 2);

  /**
   * .. function:: gsl_ran_logistic_Qinv(Q, a)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the logistic distribution
   *  with scale parameter ``a``.
   */
  ADDFUNC(gsl_cdf_logistic_Qinv, 2);

  /**
   * @file ran-pareto
   *
   * The Pareto Distribution
   * =======================
   */

  /**
   * .. function:: gsl_ran_pareto(a, b)
   *
   *  This function returns a random variate from the :index:`Pareto
   *  distribution` of order ``a``. The distribution function is,
   *
   *  .. math::
   *    p(x) dx = (a/b) / (x/b)^{a+1} dx
   *
   *  for $x \geq b$.
   */
  ADDFUNC_RANDOM(gsl_ran_pareto, 2);

  /**
   * .. function:: gsl_ran_pareto_pdf(x, a, b)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Pareto distribution with exponent ``a`` and scale ``b``, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_pareto_pdf, 3);

  /**
   * .. function:: gsl_cdf_pareto_P(x, a, b)
   */
  ADDFUNC(gsl_cdf_pareto_P, 3);

  /**
   * .. function:: gsl_cdf_pareto_Q(x, a, b)
   */
  ADDFUNC(gsl_cdf_pareto_Q, 3);

  /**
   * .. function:: gsl_cdf_pareto_Pinv(P, a, b)
   */
  ADDFUNC(gsl_cdf_pareto_Pinv, 3);

  /**
   * .. function:: gsl_cdf_pareto_Qinv(Q, a, b)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the Pareto distribution
   *  with exponent ``a`` and scale ``b``.
   */
  ADDFUNC(gsl_cdf_pareto_Qinv, 3);

  /* The spherical vector distributions are not wrapped because they return
     more than one value. */

  /**
   * @file ran-weibull
   *
   * The Weibull Distribution
   * ========================
   */

  /**
   * .. function:: gsl_ran_weibull(a, b)
   *
   *  This function returns a random variate from the :index:`Weibull
   *  distribution`. The distribution function is,
   *
   *  .. math::
   *    p(x) dx = {b \over a^b} x^{b-1} \exp(-(x/a)^b) dx
   *
   *  for $x \geq 0$.
   */
  ADDFUNC_RANDOM(gsl_ran_weibull, 2);

  /**
   * .. function:: gsl_ran_weibull_pdf(x, a, b)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Weibull distribution with scale ``a`` and exponent ``b``, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_weibull_pdf, 3);

  /**
   * .. function:: gsl_cdf_weibull_P(x, a, b)
   */
  ADDFUNC(gsl_cdf_weibull_P, 3);

  /**
   * .. function:: gsl_cdf_weibull_Q(x, a, b)
   */
  ADDFUNC(gsl_cdf_weibull_Q, 3);

  /**
   * .. function:: gsl_cdf_weibull_Pinv(P, a, b)
   */
  ADDFUNC(gsl_cdf_weibull_Pinv, 3);

  /**
   * .. function:: gsl_cdf_weibull_Qinv(Q, a, b)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the Weibull distribution
   *  with scale ``a`` and exponent ``b``.
   */
  ADDFUNC(gsl_cdf_weibull_Qinv, 3);

  /**
   * @file ran-gumbel1
   *
   * The Type-1 Gumbel Distribution
   * ==============================
   */

  /**
   * .. function:: gsl_ran_gumbel1(a, b)
   *
   *  This function returns a random variate from the :index:`Type-1 Gumbel
   *  distribution`. The Type-1 Gumbel distribution function is,
   *
   *  .. math::
   *    p(x) dx = a b \exp(-(b \exp(-ax) + ax)) dx
   *
   *  for $-\infty < x < \infty$.
   */
  ADDFUNC_RANDOM(gsl_ran_gumbel1, 2);

  /**
   * .. function:: gsl_ran_gumbel1_pdf(x, a, b)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Type-1 Gumbel distribution with parameters ``a`` and ``b``, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_gumbel1_pdf, 3);

  /**
   * .. function:: gsl_cdf_gumbel1_P(x, a, b)
   */
  ADDFUNC(gsl_cdf_gumbel1_P, 3);

  /**
   * .. function:: gsl_cdf_gumbel1_Q(x, a, b)
   */
  ADDFUNC(gsl_cdf_gumbel1_Q, 3);

  /**
   * .. function:: gsl_cdf_gumbel1_Pinv(P, a, b)
   */
  ADDFUNC(gsl_cdf_gumbel1_Pinv, 3);

  /**
   * .. function:: gsl_cdf_gumbel1_Qinv(Q, a, b)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the Type-1 Gumbel distribution
   *  with parameters ``a`` and ``b``.
   */
  ADDFUNC(gsl_cdf_gumbel1_Qinv, 3);

  /**
   * @file ran-gumbel2
   *
   * The Type-2 Gumbel Distribution
   * ==============================
   */

  /**
   * .. function:: gsl_ran_gumbel2(a, b)
   *
   *  This function returns a random variate from the :index:`Type-2 Gumbel
   *  distribution`. The Type-2 Gumbel distribution function is,
   *
   *  .. math::
   *    p(x) dx = a b x^{-a-1} \exp(-b x^{-a}) dx
   *
   *  for $-\infty < x < \infty$.
   */
  ADDFUNC_RANDOM(gsl_ran_gumbel2, 2);

  /**
   * .. function:: gsl_ran_gumbel2_pdf(x, a, b)
   *
   *  This function computes the probability density $p(x)$ at $x$ for a
   *  Type-2 Gumbel distribution with parameters ``a`` and ``b``, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_gumbel2_pdf, 3);

  /**
   * .. function:: gsl_cdf_gumbel2_P(x, a, b)
   */
  ADDFUNC(gsl_cdf_gumbel2_P, 3);

  /**
   * .. function:: gsl_cdf_gumbel2_Q(x, a, b)
   */
  ADDFUNC(gsl_cdf_gumbel2_Q, 3);

  /**
   * .. function:: gsl_cdf_gumbel2_Pinv(P, a, b)
   */
  ADDFUNC(gsl_cdf_gumbel2_Pinv, 3);

  /**
   * .. function:: gsl_cdf_gumbel2_Qinv(Q, a, b)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(x), Q(x)$ and their inverses for the Type-2 Gumbel distribution
   *  with parameters ``a`` and ``b``.
   */
  ADDFUNC(gsl_cdf_gumbel2_Qinv, 3);

  /* The Dirichlet distributions is not wrapped because it returns more
     than one value. */

  /* The general discrete distributions are not wrapped because they use
     C structures. */

  /**
   * @file ran-poisson
   *
   * The Poisson Distribution
   * ========================
   */

  /**
   * .. function:: gsl_ran_poisson(mu)
   *
   *  This function returns a random variate from the :index:`Poisson
   *  distribution` with mean ``mu``. The probability distribution
   *  for Poisson variates is,
   *
   *  .. math::
   *    p(k) = {\mu^k \over k!} \exp(-\mu)
   *
   *  for $k \geq 0$.
   */
  ADDFUNC_RANDOM(gsl_ran_poisson, 1);

  /**
   * .. function:: gsl_ran_poisson_pdf(k, mu)
   *
   *  This function computes the probability $p(k)$ of obtaining $k$ from a
   *  Poisson distribution with mean ``mu``, using the formula given above.
   */
  ADDFUNC(gsl_ran_poisson_pdf, 2);

  /**
   * .. function:: gsl_cdf_poisson_P(k, mu)
   */
  ADDFUNC(gsl_cdf_poisson_P, 2);

  /**
   * .. function:: gsl_cdf_poisson_Q(k, mu)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(k), Q(k)$ for the Poisson distribution with parameter ``mu``.
   */
  ADDFUNC(gsl_cdf_poisson_Q, 2);

  /**
   * @file ran-bernoulli
   *
   * The Bernoulli Distribution
   * ==========================
   *
   * .. index:: Bernoulli distribution
   */

  /**
   * .. function:: gsl_ran_bernoulli(p)
   *
   *  This function returns either 0 or 1, the result of a Bernoulli trial
   *  with probability ``p``. The probability distribution for a Bernoulli
   *  trial is,
   *
   *  .. math::
   *    p(0) = 1 - p \\
   *    p(1) = p
   */
  ADDFUNC_RANDOM(gsl_ran_bernoulli, 1);

  /**
   * .. function:: gsl_ran_bernoulli_pdf(k, p)
   *
   *  This function computes the probability $p(k)$ of obtaining $k$ from a
   *  Bernoulli distribution with probability parameter ``p``, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_bernoulli_pdf, 2);

  /**
   * @file ran-binomial
   *
   * The Binomial Distribution
   * =========================
   */

  /**
   * .. function:: gsl_ran_binomial(p, n)
   *
   *  This function returns a random integer from the :index:`binomial
   *  distribution`, the number of successes in ``n`` independent
   *  trials with probability ``p``. The probability distribution for
   *  binomial variates is,
   *
   *  .. math::
   *    p(k) = {n! \over k! (n-k)! } p^k (1-p)^{n-k}
   *
   *  for $0 \leq k \leq n$.
   */
  ADDFUNC_RANDOM(gsl_ran_binomial, 2);

  /**
   * .. function:: gsl_ran_binomial_pdf(k, p, n)
   *
   *  This function computes the probability $p(k)$ of obtaining $k$ from
   *  a :index:`binomial distribution` with parameters ``p`` and ``n``, using
   *  the formula given above.
   */
  ADDFUNC(gsl_ran_binomial_pdf, 3);

  /**
   * .. function:: gsl_cdf_binomial_P(k, p, n)
   */
  ADDFUNC(gsl_cdf_binomial_P, 3);

  /**
   * .. function:: gsl_cdf_binomial_Q(k, p, n)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(k), Q(k)$ for the binomial distribution with parameters
   *  ``p`` and ``n``.
   */
  ADDFUNC(gsl_cdf_binomial_Q, 3);

  /* The multinomial distributions is not wrapped because it returns more
     than one value. */

  /**
   * @file ran-negative-binomial
   *
   * The Negative Binomial Distribution
   * ==================================
   */

  /**
   * .. function:: gsl_ran_negative_binomial(p, n)
   *
   *  This function returns a random integer from the :index:`negative binomial
   *  distribution`, the number of failures occurring before ``n``
   *  successes in independent trials with probability ``p`` of success.
   *  The probability distribution for negative binomial variates is,
   *
   *  .. math::
   *    p(k) = {\Gamma(n + k) \over \Gamma(k+1) \Gamma(n) } p^n (1-p)^k
   *
   *  Note that ``n`` is not required to be an integer.
   */
  ADDFUNC_RANDOM(gsl_ran_negative_binomial, 2);

  /**
   * .. function:: gsl_ran_negative_binomial_pdf(k, p, n)
   *
   *  This function computes the probability $p(k)$ of obtaining $k$ from
   *  a negative binomial distribution with parameters ``p`` and ``n``,
   *  using the formula given above.
   */
  ADDFUNC(gsl_ran_negative_binomial_pdf, 3);

  /**
   * .. function:: gsl_cdf_negative_binomial_P(k, p, n)
   */
  ADDFUNC(gsl_cdf_negative_binomial_P, 3);

  /**
   * .. function:: gsl_cdf_negative_binomial_Q(k, p, n)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(k), Q(k)$ for the negative binomial distribution with parameters
   *  ``p`` and ``n``.
   */
  ADDFUNC(gsl_cdf_negative_binomial_Q, 3);

  /**
   * @file ran-pascal
   *
   * The Pascal Distribution
   * =======================
   */

  /**
   * .. function:: gsl_ran_pascal(p, n)
   *
   *  This function returns a random integer from the :index:`Pascal
   *  distribution`. The Pascal distribution is simply a negative
   *  binomial distribution with an integer value of ``n``.
   *
   *  .. math::
   *    p(k) = {(n + k - 1)! \over k! (n - 1)! } p^n (1-p)^k
   *
   *  for $k \geq 0$
   */
  ADDFUNC_RANDOM(gsl_ran_pascal, 2);

  /**
   * .. function:: gsl_ran_pascal_pdf(k, p, n)
   *
   *  This function computes the probability $p(k)$ of obtaining $k$ from
   *  a Pascal distribution with parameters ``p`` and ``n``, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_pascal_pdf, 3);

  /**
   * .. function:: gsl_cdf_pascal_P(k, p, n)
   */
  ADDFUNC(gsl_cdf_pascal_P, 3);

  /**
   * .. function:: gsl_cdf_pascal_Q(k, p, n)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(k), Q(k)$ for the Pascal distribution with parameters ``p``
   *  and ``n``.
   */
  ADDFUNC(gsl_cdf_pascal_Q, 3);

  /**
   * @file ran-geometric
   *
   * The Geometric Distribution
   * ==========================
   */

  /**
   * .. function:: gsl_ran_geometric(p)
   *
   *  This function returns a random integer from the :index:`geometric
   *  distribution`, the number of independent trials with
   *  probability ``p`` until the first success. The probability
   *  distribution for geometric variates is,
   *
   *  .. math::
   *    p(k) =  p (1-p)^{k-1}
   *
   *  for $k \geq 1$. Note that the distribution begins with $k=1$ with this
   *  definition. There is another convention in which the exponent $k-1$ is
   *  replaced by $k$.
   */
  ADDFUNC_RANDOM(gsl_ran_geometric, 1);

  /**
   * .. function:: gsl_ran_geometric_pdf(k, p)
   *
   *  This function computes the probability $p(k)$ of obtaining $k$ from
   *  a geometric distribution with probability parameter ``p``, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_geometric_pdf, 2);

  /**
   * .. function:: gsl_cdf_geometric_P(k, p)
   */
  ADDFUNC(gsl_cdf_geometric_P, 2);

  /**
   * .. function:: gsl_cdf_geometric_Q(k, p)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(k), Q(k)$ for the geometric distribution with parameter ``p``.
   */
  ADDFUNC(gsl_cdf_geometric_Q, 2);

  /**
   * @file ran-hypergeometric
   *
   * The Hypergeometric Distribution
   * ===============================
   */

  /**
   * .. function:: gsl_ran_hypergeometric(p, n1, n2, t)
   *
   *  This function returns a random integer from the :index:`hypergeometric
   *  distribution`. The probability distribution for hypergeometric
   *  random variates is,
   *
   *  .. math::
   *    p(k) =  C(n_1, k) C(n_2, t - k) / C(n_1 + n_2, t)
   *
   *  where $C(a,b) = a!/(b!(a-b)!)$ and $t \leq n_1 + n_2$. The domain
   *  of $k$ is $\max(0,t-n_2), ..., \min(t,n_1)$.
   *
   *  If a population contains $n_1$ elements of "type 1" and $n_2$
   *  elements of "type 2" then the hypergeometric distribution gives
   *  the probability of obtaining $k$ elements of "type 1" in $t$
   *  samples from the population without replacement.
   */
  ADDFUNC_RANDOM(gsl_ran_hypergeometric, 3);

  /**
   * .. function:: gsl_ran_hypergeometric_pdf(k, n1, n2, t)
   *
   *  This function computes the probability $p(k)$ of obtaining $k$ from
   *  a hypergeometric distribution with parameters ``n1``, ``n2``, ``t``,
   *  using the formula given above.
   */
  ADDFUNC(gsl_ran_hypergeometric_pdf, 4);

  /**
   * .. function:: gsl_cdf_hypergeometric_P(k, n1, n2, t)
   */
  ADDFUNC(gsl_cdf_hypergeometric_P, 4);

  /**
   * .. function:: gsl_cdf_hypergeometric_Q(k, n1, n2, t)
   *
   *  These functions compute the cumulative distribution functions
   *  $P(k), Q(k)$ for the hypergeometric distribution with parameters
   *  ``n1``, ``n2`` and ``t``.
   */
  ADDFUNC(gsl_cdf_hypergeometric_Q, 4);

  /**
   * @file ran-logarithmic
   *
   * The Logarithmic Distribution
   * ============================
   */

  /**
   * .. function:: gsl_ran_logarithmic(p)
   *
   *  This function returns a random integer from the :index:`logarithmic
   *  distribution`. The probability distribution for logarithmic
   *  random variates is,
   *
   *  .. math::
   *    p(k) = {-1 \over \log(1-p)} {\left(p^k \over k\right)}
   *
   *  for $k \geq 1$.
   */
  ADDFUNC_RANDOM(gsl_ran_logarithmic, 1);

  /**
   * .. function:: gsl_ran_logarithmic_pdf(k, p)
   *
   *  This function computes the probability $p(k)$ of obtaining $k$ from
   *  a logarithmic distribution with probability parameter ``p``, using the
   *  formula given above.
   */
  ADDFUNC(gsl_ran_logarithmic_pdf, 2);

  /* Shuffling and Sampling functions are not wrapped. */

  /**
   * @file ran-refs
   *
   * References and Further Reading
   * ==============================
   *
   * For an encyclopaedic coverage of the subject readers are advised to
   * consult the book *Non-Uniform Random Variate Generation* by Luc Devroye.
   * It covers every imaginable distribution and provides hundreds of
   * algorithms.
   *
   * * Luc Devroye, *Non-Uniform Random Variate Generation*, Springer-Verlag,
   *   ISBN 0-387-96305-7.
   *   Available online at http://luc.devroye.org/rnbookindex.html.
   *
   * The subject of random variate generation is also reviewed by Knuth,
   * who describes algorithms for all the major distributions.
   *
   * * Donald E. Knuth, *The Art of Computer Programming: Seminumerical
   *   Algorithms* (Vol 2, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896842.
   *
   * The Particle Data Group provides a short review of techniques for
   * generating distributions of random numbers in the "Monte Carlo"
   * section of its Annual Review of Particle Physics.
   *
   * * *Review of Particle Properties* R.M. Barnett et al., Physical Review
   *   D54, 1 (1996) http://pdg.lbl.gov/.
   *
   * The Review of Particle Physics is available online in postscript and pdf
   * format.
   *
   * An overview of methods used to compute cumulative distribution functions
   * can be found in *Statistical Computing* by W.J. Kennedy and J.E. Gentle.
   * Another general reference is *Elements of Statistical Computing* by
   * R.A. Thisted.
   *
   * * William E. Kennedy and James E. Gentle, *Statistical Computing* (1980),
   *   Marcel Dekker, ISBN 0-8247-6898-1.
   * * Ronald A. Thisted, *Elements of Statistical Computing* (1988),
   *   Chapman & Hall, ISBN 0-412-01371-1.
   *
   * The cumulative distribution functions for the Gaussian distribution
   * are based on the following papers,
   *
   * * *Rational Chebyshev Approximations Using Linear Equations*, W.J. Cody,
   *   W. Fraser, J.F. Hart. Numerische Mathematik 12, 242-251 (1968).
   * * *Rational Chebyshev Approximations for the Error Function*, W.J. Cody.
   *   Mathematics of Computation 23, n107, 631-637 (July 1969).
   */

  /**
  * @file statistics
  *
  * Statistics
  * ==========
  *
  * This chapter describes the statistical functions in the library. The basic
  * statistical functions include routines to compute the mean, variance and
  * standard deviation. More advanced functions allow you to calculate absolute
  * deviations, skewness, and kurtosis as well as the median and arbitrary
  * percentiles. The algorithms use recurrence relations to compute average
  * quantities in a stable way, without large intermediate values that might 
  * overflow.
  *
  * .. toctree::
  *    :maxdepth: 2
  *
  *    stat-mean
  *    stat-absolutedev
  *    stat-moments
  *    stat-autocorrelation
  *    stat-covariance
  *    stat-correlation
  *    stat-maxmin
  *    stat-median
  *    stat-order
  *    stat-robustlocation
  *    stat-robustscale
  *    stat-refs
  */

  /**
  * @file stat-mean
  *
  * Mean, Standard Deviation and Variance
  * =====================================
  */

  /**
  * .. function:: gsl_stats_mean(data)
  *
  * This function returns the arithmetic mean of data, a dataset
  * of length n with stride stride. The arithmetic mean, or sample mean,
  * is denoted by :math:`\hat{\mu}` and defined as,
  *
  *  .. math::
  *     \hat{\mu}= {1 \over N} \sum x_i
  *
  * where x_i are the elements of the dataset data. For samples drawn 
  * from a gaussian distribution the variance of $\Hat\mu is $\sigma^2 / $N.
  */
  ADDFUNC(gsl_stats_mean, -1);
  /**
  * .. function:: gsl_stats_variance (data)
  * This function returns the estimated, or *sample*, variance of
  * data a dataset of length *n*`. The estimated variance is denoted by 
  *:math:`\hat{\sigma^2}` and is defined by,
  *
  *  .. math::
  *    {\hat{\sigma}}^2 = {1 \over (N-1)} \sum (x_i - {\hat{\mu}})^2
  *
  * where :math:`x_i` are the elements of the dataset *data*.  Note that
  * the normalization factor of :math:`1/(N-1)` results from the derivation
  * of :math:`\hat{\sigma}^2` as an unbiased estimator of the population
  * variance :math:`\sigma^2`.  For samples drawn from a Gaussian distribution
  * the variance of :math:`\hat{\sigma}^2` itself is :math:`2 \sigma^4 / N`.
  *
  * This function computes the mean via a call to :func:`gsl_stats_mean`.  If
  * you have already computed the mean then you can pass it directly to
  * :func:`gsl_stats_variance_m`.
  */
  ADDFUNC(gsl_stats_variance, -1);

  /**
  * .. function:: gsl_stats_variance_m (data, mean)
  *
  * This function returns the sample variance of *data* relative to the
  * given value of *mean*.  The function is computed with :math:`\hat{\mu}`
  * replaced by the value of *mean* that you supply,
  * 
  *.. math:: {\hat{\sigma}}^2 = {1 \over (N-1)} \sum (x_i - mean)^2
  * 
  */
  ADDFUNC(gsl_stats_variance_m, -1);
  /**
  * 
  * .. function:: gsl_stats_sd (data)
  * .. function:: gsl_stats_sd_m(data, mean)
  * 
  * The standard deviation is defined as the square root of the variance.
  * These functions return the square root of the corresponding variance
  * functions above.
  */
  ADDFUNC(gsl_stats_sd, -1);
  ADDFUNC(gsl_stats_sd_m, -1);
  /**
  * .. function:: gsl_stats_tss (data)
  * .. function:: gsl_stats_tss_m (data, mean)
  *
  * These functions return the total sum of squares(TSS) of *data* about
  * the mean.For :func:`gsl_stats_tss_m` the user - supplied value of
  * *mean* is used, and for :func:`gsl_stats_tss` it is computed using
  * :func:`gsl_stats_mean`.
  * 
  *.. math:: {\rm TSS} = \sum(x_i - mean) ^ 2
  */
  ADDFUNC(gsl_stats_tss, -1);
  ADDFUNC(gsl_stats_tss_m, -1);
  /**
  * .. function:: gsl_stats_variance_with_fixed_mean(data, mean)
  *
  * This function computes an unbiased estimate of the variance of
  * *data* when the population mean *mean* of the underlying
  * distribution is known *a priori* .In this case the estimator for
  * the variance uses the factor :math:`1/N` and the sample mean
  * :math:`\hat{\mu}` is replaced by the known population mean :math:`\mu`,
  * 
  *.. math:: {\hat{\sigma}} ^ 2 = { 1 \over N } \sum(x_i - \mu) ^ 2
  */
  ADDFUNC(gsl_stats_variance_with_fixed_mean, -1);
  /**
  * .. function:: gsl_stats_sd_with_fixed_mean(data, mean)
  *
  * This function calculates the standard deviation of *data* for a
  * fixed population mean *mean*.  The result is the square root of the
  * corresponding variance function.
  */
  ADDFUNC(gsl_stats_sd_with_fixed_mean, -1);

  // Absolute deviation
  /**
  * @file stat-absolutedev
  *
  * Absolute deviation
  * ==================
  */
  /**
  * .. function:: gsl_stats_absdev(data)
  *
  * This function computes the absolute deviation from the mean of
  * *data*, a dataset of length *n*.  The
  * absolute deviation from the mean is defined as,
  * 
  *.. math:: absdev  = {1 \over N} \sum |x_i - {\hat{\mu}}|
  * 
  * where :math:`x_i` are the elements of the dataset *data*.  The
  * absolute deviation from the mean provides a more robust measure of the
  * width of a distribution than the variance.  This function computes the
  * mean of *data* via a call to :func:`gsl_stats_mean`.
  */
  ADDFUNC(gsl_stats_absdev, -1);
  /**
  * .. function:: gsl_stats_absdev_m(data)
  *
  * This function computes the absolute deviation of the dataset *data*
  * relative to the given value of *mean*,
  *
  *.. math:: absdev  = {1 \over N} \sum |x_i - mean|
  * 
  * This function is useful if you have already computed the mean of
  * *data* (and want to avoid recomputing it), or wish to calculate the
  * absolute deviation relative to another value (such as zero, or the
  * median).
   */
  ADDFUNC(gsl_stats_absdev_m, -1);

  // Higher moments (skewness and kurtosis)
  /**
  * @file stat-moments
  *
  * Higher moments (skewness and kurtosis)
  * ======================================
  */
  /**
  * .. function:: gsl_stats_skew(data)
  *
  * This function computes the skewness of *data*, a dataset of length
  * *n*. The skewness is defined as,
  *
  *  .. math::
  *     skew = {1 \over N} \sum 
  *       {\left( x_i - {\hat{\mu}} \over {\hat{\sigma}} \right)}^3
  * 
  * where :math:`x_i` are the elements of the dataset *data*.  The skewness
  * measures the asymmetry of the tails of a distribution.
  * 
  * The function computes the mean and estimated standard deviation of
  * *data* via calls to :func:`gsl_stats_mean` and :func:`gsl_stats_sd`.
  */
  ADDFUNC(gsl_stats_skew, -1);
  /**
  * .. function:: gsl_stats_skew_m_sd(data, mean, sd)
  *
  * This function computes the skewness of the dataset *data* using the
  * given values of the mean *mean* and standard deviation *sd*,
  *
  *  .. math::
  *     skew = {1 \over N} \sum {\left( x_i - mean \over sd \right)}^3
  *
  * where :math:`x_i` are the elements of the dataset *data*.  The skewness
  * measures the asymmetry of the tails of a distribution.
  *
  * These functions are useful if you have already computed the mean and
  * standard deviation of *data* and want to avoid recomputing them.
  */
  ADDFUNC(gsl_stats_skew_m_sd, -1);

  /**
  * .. function:: gsl_stats_kurtosis(data)
  *
  * This function computes the kurtosis of *data*. The kurtosis is defined as,
  *
  *  .. math::
  *     kurtosis = {1 \over N}  \left( \sum
  *       {\left(x_i - {\hat{\mu}} \over {\hat{\sigma}} \right)}^4 
  *       \right) - 3
  * The kurtosis measures how sharply peaked a distribution is, relative to
  * its width.  The kurtosis is normalized to zero for a Gaussian
  * distribution.
   */
  ADDFUNC(gsl_stats_kurtosis, -1);

  /**
  * .. function:: gsl_stats_kurtosis_m_sd(data, mean, sd)
  *
  * This function computes the kurtosis of the dataset *data* using the
  * given values of the mean *mean* and standard deviation *sd*,
  *
  *  .. math::
  *     kurtosis = {1 \over N}
  *         \left( \sum {\left(x_i - mean \over sd \right)}^4 \right) - 3
  * 
  * This function is useful if you have already computed the mean and
  * standard deviation of *data* and want to avoid recomputing them.
   */
  ADDFUNC(gsl_stats_kurtosis_m_sd, -1);

  // Autocorrelation
  /**
  * @file stat-autocorrelation
  *
  * Autocorrelation
  * ===============
  */
  /**
  * .. function:: gsl_stats_lag1_autocorrelation(data)
  *
  * This function computes the lag-1 autocorrelation of the dataset *data*.
  *
  *  .. math::
  *      a_1 = {\sum_{i = 2}^{n} (x_{i} - \hat{\mu}) (x_{i-1} - \hat{\mu})
  *       \over
  *       \sum_{i = 1}^{n} (x_{i} - \hat{\mu}) (x_{i} - \hat{\mu})}
  * 
  */
  ADDFUNC(gsl_stats_lag1_autocorrelation, -1);

  /**
  * .. function:: gsl_stats_lag1_autocorrelation_m (data, mean)
  *
  * This function computes the lag-1 autocorrelation of the dataset
  * *data* using the given value of the mean *mean*.
  */
  ADDFUNC(gsl_stats_lag1_autocorrelation_m, -1);


  // Covariance
  /**
  * @file stat-covariance
  *
  * Covariance
  * ==========
  */

  /**
  * .. function:: gsl_stats_covariance (data1, data2)
  *
  * This function computes the covariance of the datasets *data1* and
  * *data2* which must both be of the same length *n*.
  *
  *  .. math:: covar = {1 \over (n - 1)} \sum_{i = 1}^{n} (x_{i} - \hat{x}) (y_{i} - \hat{y})
  *
  * 
  */
  ADDFUNC(gsl_stats_covariance, -1);

  /**
  * .. function:: gsl_stats_covariance_m (data1, data2, mean1, mean2)
  *
  * This function computes the covariance of the datasets *data1* and
  * *data2* which must both be of the same length *n*.
  *
  *  .. math:: covar = {1 \over (n - 1)} \sum_{i = 1}^{n} (x_{i} - \hat{x}) (y_{i} - \hat{y})
  *
  * This function computes the covariance of the datasets *data1* and
  * *data2* using the given values of the means, *mean1* and
  * *mean2*.  This is useful if you have already computed the means of
  * *data1* and *data2* and want to avoid recomputing them.
  */
  ADDFUNC(gsl_stats_covariance_m, -1);

  // Correlation
  /**
  * @file stat-correlation
  *
  * Correlation
  * ===========
  */

  /**
  * .. function:: gsl_stats_correlation (data1, data2)
  *
  * This function efficiently computes the Pearson correlation coefficient
  * between the datasets *data1* and *data2* which must both be of
  * the same length *n*.
  *
  *  .. math:: r = {cov(x, y) \over \hat{\sigma_x} \hat{\sigma_y}} =
  *       {{1 \over n-1} \sum (x_i - \hat{x}) (y_i - \hat{y})
  *       \over
  *       \sqrt{{1 \over n-1} \sum (x_i - {\hat{x}})^2}
  *       \sqrt{{1 \over n-1} \sum (y_i - {\hat{y}})^2}
  *       }
  */ 
  ADDFUNC(gsl_stats_correlation, -1);

  /**
  * .. function:: gsl_stats_spearman (data1, data2)
  *
  *  This function computes the Spearman rank correlation coefficient between
  * the datasets *data1* and *data2* which must both be of the same
  * length *n*. The Spearman rank correlation between vectors :math:`x` and
  * :math:`y` is equivalent to the Pearson correlation between the ranked
  * vectors :math:`x_R` and :math:`y_R`, where ranks are defined to be the
  * average of the positions of an element in the ascending order of the values.
  */
  ADDFUNC(gsl_stats_spearman, -1);

  // Weighted Samples functions are not wrapped

  // Maximum and Minimum values
  /**
  * @file stat-maxmin
  *
  * Maximum and Minimum values
  * ==========================
  */

  /**
  * .. function:: gsl_stats_max (data)
  *
  * This function returns the maximum value in *data*. The maximum value is defined
  * as the value of the element :math:`x_i` which satisfies :math:`x_i \ge x_j`
  * for all :math:`j`.
  *
  * If you want instead to find the element with the largest absolute
  * magnitude you will need to apply :func:`fabs` or :func:`abs` to your data
  * before calling this function.
  */
  ADDFUNC(gsl_stats_max, -1);

  /**
  * .. function:: gsl_stats_min (data)
  *
  * This function returns the minimum value in *data*.  The minimum value is defined
  * as the value of the element :math:`x_i` which satisfies :math:`x_i \le x_j`
  * for all :math:`j`.
  *
  * If you want instead to find the element with the smallest absolute
  * magnitude you will need to apply :func:`fabs` or :func:`abs` to your data
  * before calling this function.
  */
  ADDFUNC(gsl_stats_min, -1);

  /**
  * .. function:: gsl_stats_max_index (data)
  *
  * This function returns the index of the maximum value in *data*. The maximum value is
  * defined as the value of the element :math:`x_i` which satisfies 
  * :math:`x_i \ge x_j`
  * for all :math:`j`.  When there are several equal maximum
  * elements then the first one is chosen.
  */
  ADDFUNC(gsl_stats_max_index, -1);

  /**
  * .. function:: gsl_stats_max_index (data)
  *
  * This function returns the index of the minimum value in *data*. The minimum value is
  * defined as the value of the element :math:`x_i` which satisfies
   :math:`x_i \le x_j`
  * for all :math:`j`.  When there are several equal maximum
  * elements then the first one is chosen.
  */
  ADDFUNC(gsl_stats_min_index, -1);

  // Median and Percentiles
  /**
  * @file stat-median
  *
  * Median and Percentiles
  * ======================
  * The median and percentile functions described in this section operate on
  * sorted data in :math:`O(1)` time. There is also a routine for computing
  * the median of an unsorted input array in average :math:`O(n)` time using
  * the quickselect algorithm. For convenience we use *quantiles*, measured on a scale
  * of 0 to 1, instead of percentiles (which use a scale of 0 to 100).
  */

  /**
  * .. function:: gsl_stats_median_from_sorted_data (data)
  *
  * This function returns the median value of :data:`sorted_data`.
  * The elements of the array
  * must be in ascending numerical order.  There are no checks to see
  * whether the data are sorted, so the function :func:`gsl_sort` should
  * always be used first.
  *
  * When the dataset has an odd number of elements the median is the value
  * of element :math:`(n-1)/2`.  When the dataset has an even number of
  * elements the median is the mean of the two nearest middle values,
  * elements :math:`(n-1)/2` and :math:`n/2`.  Since the algorithm for
  * computing the median involves interpolation this function always returns
  * a floating-point number, even for integer data types.
  */
  ADDFUNC(gsl_stats_median_from_sorted_data, -1);

  /**
  * .. function:: gsl_stats_median (data)
  *
  * This function returns the median value of *data*, a dataset
  * The median is found using the quickselect algorithm. The input array 
  * does not need to be sorted.
  */
  ADDFUNC(gsl_stats_median, -1);

  /**
  * .. function:: gsl_stats_quantile_from_sorted_data (data, f)
  * 
  * This function returns a quantile value of *sorted_data*. The
  * elements of the array must be in ascending numerical order.  The
  * quantile is determined by the *f*, a fraction between 0 and 1.  For
  * example, to compute the value of the 75th percentile *f* should have
  * the value 0.75.
  *
  * There are no checks to see whether the data are sorted, so the function
  * :func:`gsl_sort` should always be used first.
  *
  * The quantile is found by interpolation, using the formula
  * 
  *    .. math:: \hbox{quantile} = (1 - \delta) x_i + \delta x_{i+1}
  *
  * where :math:`i` is ``floor((n - 1)f)`` and :math:`\delta` is
  * :math:`(n-1)f - i`.
  *
  * Thus the minimum value of the array (:code:`data[1]`) is given by
  * *f* equal to zero, the maximum value (:code:`data[n]`) is
  * given by *f* equal to one and the median value is given by *f*
  * equal to 0.5.  Since the algorithm for computing quantiles involves
  * interpolation this function always returns a floating-point number, even
  * for integer data types.
  */

  ADDFUNC(gsl_stats_quantile_from_sorted_data, -1);

  // Order Statistics
  ADDFUNC(gsl_stats_select, -1);

  // Robust Location Estimates
  ADDFUNC(gsl_stats_trmean_from_sorted_data, -1);
  ADDFUNC(gsl_stats_gastwirth_from_sorted_data, -1);

  // Robust Scale  Estimates
  ADDFUNC(gsl_stats_mad0, -1);
  ADDFUNC(gsl_stats_mad, -1);
  ADDFUNC(gsl_stats_Sn0_from_sorted_data, -1);
  ADDFUNC(gsl_stats_Sn_from_sorted_data, -1);
  ADDFUNC(gsl_stats_Qn0_from_sorted_data, -1);
  ADDFUNC(gsl_stats_Qn_from_sorted_data, -1);

/**
   * @file ran-refs
   *
   * References and Further Reading
   * ==============================
   *
   * The standard reference for almost any topic in statistics is the 
   * multi-volume *Advanced Theory of Statistics* by Kendall and Stuart.
   *
   * * Maurice Kendall, Alan Stuart, and J. Keith Ord., 
   *   *The Advanced Theory of Statistics* (multiple volumes), reprinted as
   *   *Kendall’s Advanced Theory of Statistics*, Wiley, ISBN 047023380X.
   *
   * Many statistical concepts can be more easily understood by a Bayesian
   * approach. The following book by Gelman, Carlin, Stern and Rubin gives
   * a comprehensive coverage of the subject.
   *
   * * Andrew Gelman, John B. Carlin, Hal S. Stern, Donald B. Rubin. 
   *   *Bayesian Data Analysis*. Chapman & Hall, ISBN 0412039915.
   *
   * For physicists the Particle Data Group provides useful reviews of
   * Probability and Statistics in the “Mathematical Tools” section of
   * its Annual Review of Particle Physics.
   *
   * * *Review of Particle Properties* R.M. Barnett et al., Physical Review
   *   D54, 1 (1996) http://pdg.lbl.gov/.
   *
   * The Review of Particle Physics is available online in postscript and pdf
   * format.
   *
   * The following papers describe robust scale estimation,
   * * C. Croux and P. J. Rousseeuw, *Time-Efficient algorithms for two highly 
       robust estimators of scale*, Comp. Stat., Physica, Heidelberg, 1992.
   * * P. J. Rousseeuw and C. Croux, *Explicit scale estimators with high
       breakdown point*, L1-Statistical Analysis and Related Methods, pp. 77-92, 1992.
   */
}
