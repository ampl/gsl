/* movstat/apply.c
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

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>

/*
gsl_movstat_apply_accum()
  Apply moving window statistic to input vector. This is a generalized
routine to handle window endpoints and apply a given accumulator to
the input vector.

Inputs: endtype      - end point handling criteria
        x            - input vector, size n
        accum        - accumulator to apply moving window statistic
        accum_params - parameters to pass to accumulator
        y            - output vector, size n
        z            - second output vector (i.e. minmax), size n; can be NULL
        w            - workspace

Notes:
1) It is allowed to have x = y for in-place moving statistics
*/

int
gsl_movstat_apply_accum(const gsl_movstat_end_t endtype,
                        const gsl_vector * x,
                        const gsl_movstat_accum * accum,
                        void * accum_params,
                        gsl_vector * y,
                        gsl_vector * z,
                        gsl_movstat_workspace * w)
{
  if (x->size != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else if (z != NULL && x->size != z->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      const int n = (int) x->size;
      const int H = w->H; /* number of samples to left of current sample */
      const int J = w->J; /* number of samples to right of current sample */
      int i;
      double x1 = 0.0;    /* pad values for data edges */
      double xN = 0.0;
      double result[2];

      /* initialize accumulator */
      (accum->init)(w->K, w->state);

      /* pad initial window if necessary */
      if (endtype != GSL_MOVSTAT_END_TRUNCATE)
        {
          if (endtype == GSL_MOVSTAT_END_PADZERO)
            {
              x1 = 0.0;
              xN = 0.0;
            }
          else if (endtype == GSL_MOVSTAT_END_PADVALUE)
            {
              x1 = gsl_vector_get(x, 0);
              xN = gsl_vector_get(x, n - 1);
            }

          /* pad initial windows with H values */
          for (i = 0; i < H; ++i)
            (accum->insert)(x1, w->state);
        }
      else if (accum->delete_oldest == NULL) /* FIXME XXX */
        {
          /* save last K - 1 samples of x for later (needed for in-place input/output) */
          int idx1 = GSL_MAX(n - J - H, 0);
          int idx2 = n - 1;

          for (i = idx1; i <= idx2; ++i)
            w->work[i - idx1] = gsl_vector_get(x, i);
        }

      /* process input vector and fill y(0:n - J - 1) */
      for (i = 0; i < n; ++i)
        {
          double xi = gsl_vector_get(x, i);
          int idx = i - J;

          (accum->insert)(xi, w->state);

          if (idx >= 0)
            {
              (accum->get)(accum_params, result, w->state);
              gsl_vector_set(y, idx, result[0]);

              if (z != NULL)
                gsl_vector_set(z, idx, result[1]);
            }
        }

      if (endtype == GSL_MOVSTAT_END_TRUNCATE)
        {
          /* need to fill y(n-J:n-1) using shrinking windows */
          int idx1 = GSL_MAX(n - J, 0);
          int idx2 = n - 1;

          if (accum->delete_oldest == NULL)
            {
              int wsize = n - GSL_MAX(n - J - H, 0); /* size of work array */

              for (i = idx1; i <= idx2; ++i)
                {
                  int nsamp = n - GSL_MAX(i - H, 0); /* number of samples in this window */
                  int j;

                  (accum->init)(w->K, w->state);

                  for (j = wsize - nsamp; j < wsize; ++j)
                    (accum->insert)(w->work[j], w->state);

                  /* yi = acc_get [ work(i:K-2) ] */
                  (accum->get)(accum_params, result, w->state);
                  gsl_vector_set(y, i, result[0]);

                  if (z != NULL)
                    gsl_vector_set(z, i, result[1]);
                }
            }
          else
            {
              for (i = idx1; i <= idx2; ++i)
                {
                  if (i - H > 0)
                    {
                      /* delete oldest window sample as we move closer to edge */
                      (accum->delete_oldest)(w->state);
                    }

                  /* yi = acc_get [ work(i:K-2) ] */
                  (accum->get)(accum_params, result, w->state);
                  gsl_vector_set(y, i, result[0]);

                  if (z != NULL)
                    gsl_vector_set(z, i, result[1]);
                }
            }
        }
      else
        {
          /* pad final windows and fill y(n-J:n-1) */
          for (i = 0; i < J; ++i)
            {
              int idx = n - J + i;

              (accum->insert)(xN, w->state);

              if (idx >= 0)
                {
                  (accum->get)(accum_params, result, w->state);
                  gsl_vector_set(y, idx, result[0]);

                  if (z != NULL)
                    gsl_vector_set(z, idx, result[1]);
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_movstat_apply()
  Apply user-defined moving window function to input vector

Inputs: endtype - end point handling criteria
        F       - user-defined function
        x       - input vector, size n
        y       - output vector, size n
        w       - workspace
*/

int
gsl_movstat_apply(const gsl_movstat_end_t endtype, const gsl_movstat_function * F,
                  const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)
{
  int status = gsl_movstat_apply_accum(endtype, x, gsl_movstat_accum_userfunc, (void *) F, y, NULL, w);
  return status;
}
