/* bspline/greville.c
 *
 * Copyright (C) 2006, 2007, 2008, 2009 Patrick Alken
 * Copyright (C) 2008, 2011 Rhys Ulerich
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>

/* Return the location of the i-th Greville abscissa */
double
gsl_bspline_greville_abscissa (const size_t i, const gsl_bspline_workspace *w)
{
  const size_t stride = w->knots->stride;
  size_t km1 = w->spline_order - 1;
  double * data = w->knots->data + (i+1)*stride;

#if GSL_RANGE_CHECK
  if (GSL_RANGE_COND(i >= w->ncontrol))
    {
      GSL_ERROR_VAL ("Greville abscissa index out of range", GSL_EINVAL, 0);
    }
#endif

  if (km1 == 0)
    {
      /* Return interval midpoints in degenerate k = 1 case */
      km1   = 2;
      data -= stride;
    }

  return gsl_stats_mean(data, stride, km1);
}

int
gsl_bspline_init_greville (const gsl_vector *abscissae,
                           gsl_bspline_workspace *w,
                           double *abserr)
{
  /* Limited function: see https://savannah.gnu.org/bugs/index.php?34361 */

  int s;

  /* Check incoming arguments satisfy mandatory algorithmic assumptions */
  if (w->spline_order < 2)
    GSL_ERROR ("spline order must be at least 2", GSL_EINVAL);
  else if (abscissae->size < 2)
    GSL_ERROR ("abscissae->size must be at least 2", GSL_EINVAL);
  else if (w->nbreak != abscissae->size - w->spline_order + 2)
    GSL_ERROR ("w->nbreak must equal abscissae->size - spline_order + 2", GSL_EINVAL);

  if (w->nbreak == 2)
    {
      /* No flexibility in abscissae values possible in this degenerate case */
      s = gsl_bspline_init_uniform (
              gsl_vector_get (abscissae, 0),
              gsl_vector_get (abscissae, abscissae->size - 1), w);
    }
  else
    {
      double * storage;
      gsl_matrix_view A;
      gsl_vector_view tau, b, x, r;
      size_t i, j;

      /* Constants derived from the B-spline workspace and abscissae details */
      const size_t km2    = w->spline_order - 2;
      const size_t M      = abscissae->size - 2;
      const size_t N      = w->nbreak - 2;
      const double invkm1 = 1.0 / (w->spline_order - 1);

      /* Allocate working storage and prepare multiple, zero-filled views */
      storage = (double *) calloc (M*N + 2*N + 2*M, sizeof (double));
      if (storage == 0)
        GSL_ERROR ("failed to allocate working storage", GSL_ENOMEM);
      A   = gsl_matrix_view_array (storage, M, N);
      tau = gsl_vector_view_array (storage + M*N,             N);
      b   = gsl_vector_view_array (storage + M*N + N,         M);
      x   = gsl_vector_view_array (storage + M*N + N + M,     N);
      r   = gsl_vector_view_array (storage + M*N + N + M + N, M);

      /* Build matrix from interior breakpoints to interior Greville abscissae.
       * For example, when spline_order = 4 and nbreak = 7 the matrix is
       *   [   1,      0,      0,      0,      0;
       *     2/3,    1/3,      0,      0,      0;
       *     1/3,    1/3,    1/3,      0,      0;
       *       0,    1/3,    1/3,    1/3,      0;
       *       0,      0,    1/3,    1/3,    1/3;
       *       0,      0,      0,    1/3,    2/3;
       *       0,      0,      0,      0,      1  ]
       * but only center formed as first/last breakpoint is known.
       */
      for (j = 0; j < N; ++j)
        for (i = 0; i <= km2; ++i)
          gsl_matrix_set (&A.matrix, i+j, j, invkm1);

      /* Copy interior collocation points from abscissae into b */
      for (i = 0; i < M; ++i)
        gsl_vector_set (&b.vector, i, gsl_vector_get (abscissae, i+1));

      /* Adjust b to account for constraint columns not stored in A */
      for (i = 0; i < km2; ++i)
        {
          double * const v = gsl_vector_ptr (&b.vector, i);
          *v -= (1 - (i+1)*invkm1) * gsl_vector_get (abscissae, 0);
        }
      for (i = 0; i < km2; ++i)
        {
          double * const v = gsl_vector_ptr (&b.vector, M - km2 + i);
          *v -= (i+1)*invkm1 * gsl_vector_get (abscissae, abscissae->size - 1);
        }

      /* Perform linear least squares to determine interior breakpoints */
      s =  gsl_linalg_QR_decomp (&A.matrix, &tau.vector)
        || gsl_linalg_QR_lssolve (&A.matrix, &tau.vector,
                                  &b.vector, &x.vector, &r.vector);
      if (s)
        {
          free (storage);
          return s;
        }

      /* "Expand" solution x by adding known first and last breakpoints. */
      x = gsl_vector_view_array_with_stride (
          gsl_vector_ptr (&x.vector, 0) - x.vector.stride,
          x.vector.stride, x.vector.size + 2);
      gsl_vector_set (&x.vector, 0, gsl_vector_get (abscissae, 0));
      gsl_vector_set (&x.vector, x.vector.size - 1,
                      gsl_vector_get (abscissae, abscissae->size - 1));

      /* Finally, initialize workspace knots using the now-known breakpoints */
      s = gsl_bspline_init_augment (&x.vector, w);
      free (storage);
    }

  /* Sum absolute errors in the resulting vs requested interior abscissae */
  /* Provided as a fit quality metric which may be monitored by callers */
  if (!s && abserr)
    {
      size_t i;
      *abserr = 0;
      for (i = 1; i < abscissae->size - 1; ++i)
        *abserr += fabs (   gsl_bspline_greville_abscissa (i, w)
                          - gsl_vector_get (abscissae, i) );
    }

  return s;
}
