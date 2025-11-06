/* bspline/bspline.c
 *
 * Copyright (C) 2006-2020 Patrick Alken
 * Copyright (C) 2008 Rhys Ulerich
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
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_statistics.h>

/*
 * This module contains routines related to calculating B-splines.
 * The algorithms used are described in
 *
 * [1] Carl de Boor, "A Practical Guide to Splines", Springer
 *     Verlag, 1978.
 *
 * The bspline_pppack_* internal routines contain code adapted from
 *
 * [2] "PPPACK - Piecewise Polynomial Package",
 *     http://www.netlib.org/pppack/
 */

/*
gsl_bspline_alloc()
  Allocate space for a bspline workspace. The size of the
workspace is O(5*order + nbreak)

Inputs: spline_order - spline order (cubic = 4)
        nbreak       - number of breakpoints

Return: pointer to workspace
*/

gsl_bspline_workspace *
gsl_bspline_alloc (const size_t spline_order, const size_t nbreak)
{
  gsl_bspline_workspace *w;

  if (spline_order == 0)
    {
      GSL_ERROR_NULL ("spline order must be at least 1", GSL_EINVAL);
    }

  if (nbreak < 2)
    {
      GSL_ERROR_NULL ("nbreak must be at least 2", GSL_EINVAL);
    }

  w = calloc (1, sizeof (gsl_bspline_workspace));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace",
        GSL_ENOMEM);
    }

  w->spline_order = spline_order;
  w->nbreak = nbreak;
  w->ncontrol = w->nbreak + spline_order - 2;
  w->icache = 0;

  w->knots = gsl_vector_alloc (w->ncontrol + spline_order);
  if (w->knots == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL ("failed to allocate space for knots vector",
                      GSL_ENOMEM);
    }

  w->deltal = gsl_vector_alloc (spline_order);
  if (w->deltal == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL ("failed to allocate space for deltal vector",
                      GSL_ENOMEM);
    }

  w->deltar = gsl_vector_alloc (spline_order);
  if (w->deltar == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL ("failed to allocate space for deltar vector",
                      GSL_ENOMEM);
    }

  w->B = gsl_vector_alloc (spline_order);
  if (w->B == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL ("failed to allocate space for temporary spline vector",
                      GSL_ENOMEM);
    }

  w->XTX = gsl_matrix_alloc (w->ncontrol, spline_order);
  if (w->XTX == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL("failed to allocate space for XTX matrix", GSL_ENOMEM);
    }

  w->R = gsl_matrix_alloc (w->ncontrol, w->ncontrol);
  if (w->R == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL("failed to allocate space for R matrix", GSL_ENOMEM);
    }

  w->work = gsl_vector_alloc(3 * w->ncontrol);
  if (w->work == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL("failed to allocate space for work array", GSL_ENOMEM);
    }

  w->A = gsl_matrix_alloc (spline_order, spline_order);
  if (w->A == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL("failed to allocate space for derivative work matrix",
                     GSL_ENOMEM);
    }

  w->dB = gsl_matrix_alloc (spline_order, 2 * (spline_order + 1));
  if (w->dB == 0)
    {
      gsl_bspline_free (w);
      GSL_ERROR_NULL ("failed to allocate space for temporary derivative matrix",
                      GSL_ENOMEM);
    }

  return w;
}

gsl_bspline_workspace *
gsl_bspline_alloc_ncontrol (const size_t spline_order, const size_t ncontrol)
{
  if (ncontrol < spline_order)
    {
      GSL_ERROR_NULL("ncontrol must be at least spline_order", GSL_EINVAL);
    }
  else
    {
      size_t nbreak = ncontrol - spline_order + 2;
      return gsl_bspline_alloc(spline_order, nbreak);
    }
}

/*
gsl_bspline_free()
  Free a gsl_bspline_workspace.

Inputs: w - workspace to free

Return: none
*/

void
gsl_bspline_free (gsl_bspline_workspace * w)
{
  if (w->knots)
    gsl_vector_free (w->knots);

  if (w->deltal)
    gsl_vector_free (w->deltal);

  if (w->deltar)
    gsl_vector_free (w->deltar);

  if (w->B)
    gsl_vector_free (w->B);

  if (w->XTX)
    gsl_matrix_free(w->XTX);

  if (w->R)
    gsl_matrix_free(w->R);

  if (w->work)
    gsl_vector_free (w->work);

  if (w->A)
    gsl_matrix_free(w->A);

  if (w->dB)
    gsl_matrix_free(w->dB);

  free (w);
}


/* return number of B-spline control points */
size_t
gsl_bspline_ncontrol (const gsl_bspline_workspace * w)
{
  return w->ncontrol;
}

/* return spline order */
size_t
gsl_bspline_order (const gsl_bspline_workspace * w)
{
  return w->spline_order;
}

/* Return number of breakpoints */
size_t
gsl_bspline_nbreak (const gsl_bspline_workspace * w)
{
  return w->nbreak;
}

/* Return the location of the i-th breakpoint */
double
gsl_bspline_breakpoint (const size_t i, const gsl_bspline_workspace * w)
{
  if (i >= w->nbreak)
    {
      GSL_ERROR_VAL("invalid breakpoint index", GSL_EDOM, 0.0);
    }
  else
    {
      size_t j = i + w->spline_order - 1;
      return gsl_vector_get (w->knots, j);
    }
}

/*
gsl_bspline_init_augment()
  Initialize knots vector using given interior knots (breakpoints):

   knots(1:k) = tau(1)
   knots(k+1:k+l-1) = tau(i), i = 2 .. l
   knots(n+1:n+k) = tau(l + 1)

where l is the number of polynomial pieces (l = nbreak - 1) and
   n = k + l - 1
(using matlab syntax for the arrays)

The repeated knots at the beginning and end of the interval
correspond to the continuity condition there. See pg. 119
of [1].

Inputs: tau - interior knots, length n - k + 2
        w   - bspline workspace

Return: success or error
*/

int
gsl_bspline_init_augment (const gsl_vector * tau, gsl_bspline_workspace * w)
{
  if (tau->size != w->nbreak)
    {
      GSL_ERROR ("tau vector size does not match workspace", GSL_EBADLEN);
    }
  else
    {
      const double b1 = gsl_vector_get(tau, 0);
      const double bn = gsl_vector_get(tau, w->nbreak - 1);
      size_t i;

      for (i = 0; i < w->spline_order; i++)
        gsl_vector_set (w->knots, i, b1);

      for (i = 1; i < w->nbreak - 1; i++)
        {
          gsl_vector_set (w->knots, w->spline_order - 1 + i,
                          gsl_vector_get (tau, i));
        }

      for (i = w->ncontrol; i < w->ncontrol + w->spline_order; i++)
        gsl_vector_set (w->knots, i, bn);

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_init_uniform()
  Initialize knots vector using uniformly spaced breakpoints on the interval [a,b]
'a' is the position of the first breakpoint and 'b' is the position of the last
breakpoint.

Inputs: a - left side of interval
        b - right side of interval
        w - bspline workspace

Return: success or error

Notes: 1) w->knots is modified to contain the uniformly spaced
          knots

       2) The knots vector is set up as follows (using octave syntax):

          knots(1:k) = a
          knots(k+1:k+l-1) = a + i*delta, i = 1 .. l - 1
          knots(n+1:n+k) = b
*/

int
gsl_bspline_init_uniform (const double a, const double b,
                          gsl_bspline_workspace * w)
{
  size_t i;     /* looping */
  double delta; /* interval spacing */
  double x;

  delta = (b - a) / (w->nbreak - 1.0);

  for (i = 0; i < w->spline_order; i++)
    gsl_vector_set (w->knots, i, a);

  x = a + delta;
  for (i = 0; i < w->nbreak - 2; i++)
    {
      gsl_vector_set (w->knots, w->spline_order + i, x);
      x += delta;
    }

  for (i = w->ncontrol; i < w->ncontrol + w->spline_order; i++)
    gsl_vector_set (w->knots, i, b);

  return GSL_SUCCESS;
}

/*
gsl_bspline_init_periodic()
  Initialize uniform periodic knots vector on the interval [a,b]

Inputs: a - left side of interval
        b - right side of interval
        w - bspline workspace

Return: success or error

Notes:
1) The resulting B-spline is assumed periodic on [a,b]

2) w->knots is modified to contain the uniformly spaced periodic knots

3) The knots vector is set up as follows (using octave syntax):

t_i = (i-k)/(n-k+1) (b-a) + a,    i = 1, ..., n + k
*/

int
gsl_bspline_init_periodic (const double a, const double b,
                           gsl_bspline_workspace * w)
{
  const double spline_order = (double) w->spline_order;
  const double denom = (double) w->ncontrol - spline_order + 1.0;
  size_t i;

  for (i = 0; i < w->knots->size; ++i)
    {
      double di = (double) i;
      double ui = (di - spline_order + 1.0) / denom;
      double ti = (b - a) * ui + a;

      gsl_vector_set(w->knots, i, ti);
    }

  return GSL_SUCCESS;
}

int
gsl_bspline_init (const gsl_vector * t, gsl_bspline_workspace * w)
{
  if (t->size != w->knots->size)
    {
      GSL_ERROR ("t vector size does not match workspace", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_memcpy(w->knots, t);
      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_proj_rhs()
  Compute the rhs vector for a projection of a function f(x) onto
the B-spline basis

f(x) =~ proj(f) = \sum_i c_i B_i(x)

The c_i can be found by solving:

G c = y

where y_i = \int f(x) B_i(x) dx

Inputs: F - function to be projected
        y - (output) rhs vector y, length ncontrol
        w - workspace

Notes:
1) The algorithm uses Gauss-Legendre quadrature on each knot interval

2) If f(x) is a piecewise polynomial of order k-1 or less, the projection
will be exact, since f lies in the space spanned by the B_i(x)
*/

int
gsl_bspline_proj_rhs(const gsl_function * F, gsl_vector * y, gsl_bspline_workspace * w)
{
  if (y->size != w->ncontrol)
    {
      GSL_ERROR ("size of y vector must match number of control points", GSL_EBADLEN);
    }
  else
    {
      const size_t k = gsl_bspline_order(w);
      const size_t m = (k + 1) >> 1;
      gsl_integration_glfixed_table * gltable = gsl_integration_glfixed_table_alloc(k);
      gsl_vector_view N1 = gsl_matrix_column(w->dB, 0);
      gsl_vector_view N2 = gsl_matrix_column(w->dB, 1);
      size_t ustart = 0;
      size_t nu, istart;
      size_t i;

      gsl_vector_set_zero(y);

      if (k & 1)
        {
          ustart = 1;

          for (nu = 0; nu < w->knots->size - 1; ++nu)
            {
              double t0 = gsl_vector_get(w->knots, nu);
              double t1 = gsl_vector_get(w->knots, nu + 1);

              if (t1 > t0)
                {
                  double A = 0.5 * (t1 - t0);
                  double B = t0 + A; /* 0.5 * (t1 + t0) */
                  double w0 = A * gltable->w[0];
                  double f0 = GSL_FN_EVAL(F, B);

                  gsl_bspline_basis(B, &N1.vector, &istart, w);

                  for (i = 0; i < k; ++i)
                    {
                      double Ni = gsl_vector_get(&N1.vector, i);
                      double * yi = gsl_vector_ptr(y, istart + i);
                      *yi += w0 * f0 * Ni;
                    }
                }
            }
        }

      for (nu = 0; nu < w->knots->size - 1; ++nu)
        {
          double t0 = gsl_vector_get(w->knots, nu);
          double t1 = gsl_vector_get(w->knots, nu + 1);

          if (t1 > t0)
            {
              double A = 0.5 * (t1 - t0);
              double B = t0 + A; /* 0.5 * (t1 + t0) */
              size_t u;

              for (u = ustart; u < m; ++u)
                {
                  double Ax = A * gltable->x[u];
                  double wu = A * gltable->w[u];
                  double fp, fm;

                  gsl_bspline_basis(B + Ax, &N1.vector, &istart, w);
                  gsl_bspline_basis(B - Ax, &N2.vector, &istart, w);

                  fp = GSL_FN_EVAL(F, B + Ax);
                  fm = GSL_FN_EVAL(F, B - Ax);

                  for (i = 0; i < k; ++i)
                    {
                      double N1i = gsl_vector_get(&N1.vector, i);
                      double N2i = gsl_vector_get(&N2.vector, i);
                      double * yi = gsl_vector_ptr(y, istart + i);

                      *yi += wu * (fp * N1i + fm * N2i);
                    }
                }
            }
        }

      gsl_integration_glfixed_table_free(gltable);

      return GSL_SUCCESS;
    }
}
