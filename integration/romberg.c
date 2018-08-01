/* integration/romberg.c
 * 
 * Copyright (C) 2018 Patrick Alken
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

/* the code in this module performs Romberg integration */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#define ROMBERG_PRINT_ROW(i, r) \
  do { \
    size_t jj; \
    fprintf(stderr, "R[%zu] = ", i); \
    for (jj = 0; jj <= i; ++jj) \
      fprintf(stderr, "%.8e ", r[jj]); \
    fprintf(stderr, "\n"); \
  } while (0)

gsl_integration_romberg_workspace *
gsl_integration_romberg_alloc(const size_t n)
{
  gsl_integration_romberg_workspace *w;

  /* check inputs */
  if (n < 1)
    {
      GSL_ERROR_VAL ("workspace size n must be at least 1", GSL_EDOM, 0);
    }

  w = calloc(1, sizeof(gsl_integration_romberg_workspace));
  if (w == NULL)
    {
      GSL_ERROR_VAL ("unable to allocate workspace", GSL_ENOMEM, 0);
    }

  /* ceiling on n, since the number of points is 2^n + 1 */
  w->n = GSL_MIN(n, 30);

  w->work1 = malloc(w->n * sizeof(double));
  if (w->work1 == NULL)
    {
      gsl_integration_romberg_free(w);
      GSL_ERROR_VAL ("unable to allocate previous row", GSL_ENOMEM, 0);
    }

  w->work2 = malloc(w->n * sizeof(double));
  if (w->work2 == NULL)
    {
      gsl_integration_romberg_free(w);
      GSL_ERROR_VAL ("unable to allocate current row", GSL_ENOMEM, 0);
    }

  return w;
}

void
gsl_integration_romberg_free(gsl_integration_romberg_workspace * w)
{
  if (w->work1)
    free(w->work1);

  if (w->work2)
    free(w->work2);

  free(w);
}

int
gsl_integration_romberg(const gsl_function * f, const double a, const double b,
                        const double epsabs, const double epsrel, double * result,
                        size_t * neval, gsl_integration_romberg_workspace * w)
{
  if (epsabs < 0.0)
    {
      GSL_ERROR("epsabs must be non-negative", GSL_EDOM);
    }
  else if (epsrel < 0.0)
    {
      GSL_ERROR("epsrel must be non-negative", GSL_EDOM);
    }
  else
    {
      const size_t n = w->n;
      double *Rp = &(w->work1[0]); /* previous row */
      double *Rc = &(w->work2[0]); /* current row */
      double *Rtmp;
      double h = 0.5 * (b - a);    /* step size */
      size_t i;

      /* R(0,0) */
      Rp[0] = h * (GSL_FN_EVAL(f, a) + GSL_FN_EVAL(f, b));
      *neval = 2;

      /*ROMBERG_PRINT_ROW((size_t) 0, Rp);*/

      for (i = 1; i < n; ++i)
        {
          size_t j;
          double sum = 0.0;
          double err;
          double four_j;         /* 4^j */
          size_t two_i = 1 << i; /* 2^i */

          for (j = 1; j < two_i; j += 2)
            {
              sum += GSL_FN_EVAL(f, a + j * h);
              ++(*neval);
            }

          /* R(i,0) */
          Rc[0] = sum * h + 0.5 * Rp[0];

          four_j = 4.0;
          for (j = 1; j <= i; ++j)
            {
              Rc[j] = (four_j * Rc[j - 1] - Rp[j - 1]) / (four_j - 1.0);
              four_j *= 4.0;
            }

          /*ROMBERG_PRINT_ROW(i, Rc);*/

          /*
           * compute difference between current and previous result and
           * check for convergence
           */
          err = fabs(Rc[i] - Rp[i - 1]);
          if ((err < epsabs) || (err < epsrel * fabs(Rc[i])))
            {
              *result = Rc[i];
              return GSL_SUCCESS;
            }

          /* swap Rp and Rc */
          Rtmp = Rp;
          Rp = Rc;
          Rc = Rtmp;

          h *= 0.5;
        }

      /* did not converge - return best guess */
      *result = Rp[n - 1];

      return GSL_EMAXITER;
    }
}
