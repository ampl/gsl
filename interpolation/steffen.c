/* interpolation/steffen.c
 * 
 * Copyright (C) 2014 Jean-François Caron
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  J.-F. Caron
 *
 * This interpolation method is taken from 
 * M.Steffen, "A simple method for monotonic interpolation in one dimension",
 * Astron. Astrophys. 239, 443-450 (1990).
 *
 * This interpolation method guarantees monotonic interpolation functions between
 * the given data points.  A consequence of this is that extremal values can only
 * occur at the data points.  The interpolating function and its first derivative
 * are guaranteed to be continuous, but the second derivative is not.
 *
 * The implementation is modelled on the existing Akima interpolation method
 * previously included in GSL by Gerard Jungman.
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "integ_eval.h"
#include <gsl/gsl_interp.h>

typedef struct
{
  double * a;       /* eqs 2-5 of paper */
  double * b;
  double * c;
  double * d;

  double * y_prime; /* eq 11 of paper */
} steffen_state_t;

static void steffen_free (void * vstate);
static double steffen_copysign(const double x, const double y);

static void *
steffen_alloc (size_t size)
{
  steffen_state_t *state;
  
  state = (steffen_state_t *) calloc (1, sizeof (steffen_state_t));
  
  if (state == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for state", GSL_ENOMEM);
    }

  state->a = (double *) malloc (size * sizeof (double));
  
  if (state->a == NULL)
    {
      steffen_free(state);
      GSL_ERROR_NULL("failed to allocate space for a", GSL_ENOMEM);
    }
  
  state->b = (double *) malloc (size * sizeof (double));
  
  if (state->b == NULL)
    {
      steffen_free(state);
      GSL_ERROR_NULL("failed to allocate space for b", GSL_ENOMEM);
    }
  
  state->c = (double *) malloc (size * sizeof (double));
  
  if (state->c == NULL)
    {
      steffen_free(state);
      GSL_ERROR_NULL("failed to allocate space for c", GSL_ENOMEM);
    }
  
  state->d = (double *) malloc (size * sizeof (double));
  
  if (state->d == NULL)
    {
      steffen_free(state);
      GSL_ERROR_NULL("failed to allocate space for d", GSL_ENOMEM);
    }

  state->y_prime = (double *) malloc (size * sizeof (double));
  if (state->y_prime == NULL)
    {
      steffen_free(state);
      GSL_ERROR_NULL("failed to allocate space for y_prime", GSL_ENOMEM);
    }

  return state;
}

static int
steffen_init (void * vstate, const double x_array[],
              const double y_array[], size_t size)
{
  steffen_state_t *state = (steffen_state_t *) vstate;
  size_t i;
  double *a = state->a;
  double *b = state->b;
  double *c = state->c;
  double *d = state->d;
  double *y_prime = state->y_prime;

  /*
   * first assign the interval and slopes for the left boundary.
   * We use the "simplest possibility" method described in the paper
   * in section 2.2
   */
  double h0 = (x_array[1] - x_array[0]);
  double s0 = (y_array[1] - y_array[0]) / h0;

  y_prime[0] = s0;

  /* Now we calculate all the necessary s, h, p, and y' variables 
     from 1 to N-2 (0 to size - 2 inclusive) */
  for (i = 1; i < (size - 1); i++)
    {
      double pi;

      /* equation 6 in the paper */
      double hi = (x_array[i+1] - x_array[i]);
      double him1 = (x_array[i] - x_array[i - 1]);

      /* equation 7 in the paper */
      double si = (y_array[i+1] - y_array[i]) / hi;
      double sim1 = (y_array[i] - y_array[i - 1]) / him1;

      /* equation 8 in the paper */
      pi = (sim1*hi + si*him1) / (him1 + hi);

      /* This is a C equivalent of the FORTRAN statement below eqn 11 */
      y_prime[i] = (steffen_copysign(1.0,sim1) + steffen_copysign(1.0,si)) *
                    GSL_MIN(fabs(sim1),
                            GSL_MIN(fabs(si), 0.5*fabs(pi))); 
    }

  /*
   * we also need y' for the rightmost boundary; we use the
   * "simplest possibility" method described in the paper in
   * section 2.2
   *
   * y' = s_{n-1}
   */
  y_prime[size-1] = (y_array[size - 1] - y_array[size - 2]) /
                    (x_array[size - 1] - x_array[size - 2]);

  /* Now we can calculate all the coefficients for the whole range. */
  for (i = 0; i < (size - 1); i++)
    {
      double hi = (x_array[i+1] - x_array[i]);
      double si = (y_array[i+1] - y_array[i]) / hi;

      /* These are from equations 2-5 in the paper. */
      a[i] = (y_prime[i] + y_prime[i+1] - 2*si) / hi / hi;
      b[i] = (3*si - 2*y_prime[i] - y_prime[i+1]) / hi;
      c[i] = y_prime[i];
      d[i] = y_array[i];
    }

  return GSL_SUCCESS;
}

static void
steffen_free (void * vstate)
{
  steffen_state_t *state = (steffen_state_t *) vstate;

  RETURN_IF_NULL(state);

  if (state->a)
    free (state->a);

  if (state->b)
    free (state->b);

  if (state->c)
    free (state->c);

  if (state->d)
    free (state->d);

  if (state->y_prime)
    free (state->y_prime);

  free (state);
}

static int
steffen_eval (const void * vstate,
              const double x_array[], const double y_array[], size_t size,
              double x, gsl_interp_accel * a, double *y)
{
  const steffen_state_t *state = (const steffen_state_t *) vstate;

  size_t index;
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  {
    const double x_lo = x_array[index];
    const double delx = x - x_lo;
    const double a = state->a[index];
    const double b = state->b[index];
    const double c = state->c[index];
    const double d = state->d[index];
    /* Use Horner's scheme for efficient evaluation of polynomials. */
    /* *y = a*delx*delx*delx + b*delx*delx + c*delx + d; */
    *y = d + delx*(c + delx*(b + delx*a));

    return GSL_SUCCESS;
  }
}

static int
steffen_eval_deriv (const void * vstate,
                    const double x_array[], const double y_array[], size_t size,
                    double x, gsl_interp_accel * a, double *dydx)
{
  const steffen_state_t *state = (const steffen_state_t *) vstate;

  size_t index;

  /* DISCARD_POINTER(y_array); /\* prevent warning about unused parameter *\/ */
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  {
    double x_lo = x_array[index];
    double delx = x - x_lo;
    double a = state->a[index];
    double b = state->b[index];
    double c = state->c[index];
    /*double d = state->d[index];*/
    /* *dydx = 3*a*delx*delx*delx + 2*b*delx + c; */
    *dydx = c + delx*(2*b + delx*3*a);
    return GSL_SUCCESS;
  }
}

static int
steffen_eval_deriv2 (const void * vstate,
                     const double x_array[], const double y_array[], size_t size,
                     double x, gsl_interp_accel * a, double *y_pp)
{
  const steffen_state_t *state = (const steffen_state_t *) vstate;

  size_t index;

  /* DISCARD_POINTER(y_array); /\* prevent warning about unused parameter *\/ */

  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  {
    const double x_lo = x_array[index];
    const double delx = x - x_lo;
    const double a = state->a[index];
    const double b = state->b[index];
    *y_pp = 6*a*delx + 2*b;
    return GSL_SUCCESS;
  }
}

static int
steffen_eval_integ (const void * vstate,
                    const double x_array[], const double y_array[], size_t size,
                    gsl_interp_accel * acc, double a, double b,
                    double * result)
{
  /* a and b are the boundaries of the integration. */
  
  const steffen_state_t *state = (const steffen_state_t *) vstate;

  size_t i, index_a, index_b;

  /* Find the data points in the x_array that are nearest to the desired */
  /* a and b integration boundaries. */

  if (acc != 0)
    {
      index_a = gsl_interp_accel_find (acc, x_array, size, a);
      index_b = gsl_interp_accel_find (acc, x_array, size, b);
    }
  else
    {
      index_a = gsl_interp_bsearch (x_array, a, 0, size - 1);
      index_b = gsl_interp_bsearch (x_array, b, 0, size - 1);
    }
  
  *result = 0.0;

  /* Iterate over all the segments between data points and sum the */
  /* contributions into result. */
  for(i=index_a; i<=index_b; i++) 
    {
      const double x_hi = x_array[i + 1];
      const double x_lo = x_array[i];
      const double dx = x_hi - x_lo;
      if(dx != 0.0) 
        {
          /*
           * check if we are at a boundary point, so take the
           * a and b parameters instead of the data points.
           */
          double x1 = (i == index_a) ? a-x_lo : 0.0;
          double x2 = (i == index_b) ? b-x_lo : x_hi-x_lo;

          *result += (1.0/4.0)*state->a[i]*(x2*x2*x2*x2 - x1*x1*x1*x1)
                    +(1.0/3.0)*state->b[i]*(x2*x2*x2 - x1*x1*x1)
                    +(1.0/2.0)*state->c[i]*(x2*x2 - x1*x1)
                    +state->d[i]*(x2-x1);
        }
      else /* if the interval was zero, i.e. consecutive x values in data. */
        {
          *result = 0.0;
          return GSL_EINVAL;
        }
    }
  
  return GSL_SUCCESS;
}

static double
steffen_copysign(const double x, const double y)
{
  if ((x < 0 && y > 0) || (x > 0 && y < 0))
    return -x;

  return x;
}

static const gsl_interp_type steffen_type = 
{
  "steffen", 
  3,
  &steffen_alloc,
  &steffen_init,
  &steffen_eval,
  &steffen_eval_deriv,
  &steffen_eval_deriv2,
  &steffen_eval_integ,
  &steffen_free
};

const gsl_interp_type * gsl_interp_steffen = &steffen_type;
