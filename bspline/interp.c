/* bspline/interp.c
 *
 * Copyright (C) 2022 Patrick Alken
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

/* routines related to B-spline interpolation */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>

/*
gsl_bspline_init_interp()
  Initialize a knot vector suitable for interpolating at abscissa values x

Inputs: x - interpolation abscissa values, nondecreasing, length ncontrol
            condition: x(i) < x(i+k-1)
        w - workspace

Return: success/error

Notes:
1) w->knots will satisfy the Schoenberg-Whitney conditions, and is initialized as follows:

knots[0:k-1]   = x[0]
knots[n:n+k-1] = x[n-1]
*/

int
gsl_bspline_init_interp (const gsl_vector * x, gsl_bspline_workspace * w)
{
  const size_t n = x->size;

  if (n != w->ncontrol)
    {
      GSL_ERROR("x length must be equal to ncontrol", GSL_EBADLEN);
    }
  else
    {
      const double x0 = gsl_vector_get(x, 0);
      const double x1 = gsl_vector_get(x, n - 1);
      const size_t k = w->spline_order;
      size_t i;

      /* set first and last knots to x(1) and x(n) with multiplicity k */
      for (i = 0; i < k; ++i)
        {
          gsl_vector_set(w->knots, i, x0);
          gsl_vector_set(w->knots, n + k - i - 1, x1);
        }

      if (k == 1)
        {
          /* use midpoints between abscissa values */
          for (i = 0; i < n - 1; ++i)
            {
              double xi = gsl_vector_get(x, i);
              double xip1 = gsl_vector_get(x, i + 1);

              if (xip1 < xi)
                {
                  GSL_ERROR("x must be non-decreasing", GSL_EINVAL);
                }

              gsl_vector_set(w->knots, i + 1, 0.5 * (xip1 + xi));
            }
        }
      else
        {
          for (i = 0; i < n - k; ++i)
            {
              gsl_vector_const_view v = gsl_vector_const_subvector(x, i + 1, k - 1);
              double mean = gsl_stats_mean(v.vector.data, v.vector.stride, k - 1);
              gsl_vector_set(w->knots, i + k, mean);
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_init_hermite()
  Initialize a knot vector for Hermite interpolation:

Given a vector of interpolation sites,

x = [ x_1, x_2, ..., x_n ]

the knot vector is:

t = [ x_0, x_0, x_1, ..., x_1, x_2, ..., x_n, ..., x_n, x_{n+1}, ..., x_{n+1} ]

where:

x_0     = x_1
x_{n+1} = x_n

and each x_i appears (m+1) times, where m is the order of derivatives we want to
match (m = 1 for the usual Hermite interpolation)

Inputs: nderiv - derivative order (m = 1 for usual Hermite interpolation)
                 (note: it is required that spline_order = 2*nderiv + 2)
        x      - vector of interpolation sites
        w      - workspace

Notes:
1) On output, w->knots is initialized

2) It is required that spline_order = 2*nderiv + 2
*/

int
gsl_bspline_init_hermite(const size_t nderiv, const gsl_vector * x, gsl_bspline_workspace * w)
{
  const size_t n = x->size;
  const size_t nknots = w->knots->size;

  if (w->spline_order != 2*nderiv + 2)
    {
      GSL_ERROR ("spline_order must be 2*nderiv+2", GSL_EDOM);
    }
  else if (nknots != (nderiv + 1) * (n + 2))
    {
      GSL_ERROR ("nknots must equal (nderiv+1)(n+2)", GSL_EBADLEN);
    }
  else
    {
      double x_i;
      size_t i, j;
      size_t idx = 0;

      x_i = gsl_vector_get(x, 0);
      for (j = 0; j <= nderiv; ++j)
        gsl_vector_set(w->knots, idx++, x_i);

      for (i = 0; i < n; ++i)
        {
          x_i = gsl_vector_get(x, i);

          for (j = 0; j <= nderiv; ++j)
            gsl_vector_set(w->knots, idx++, x_i);
        }

      x_i = gsl_vector_get(x, n - 1);
      for (j = 0; j <= nderiv; ++j)
        gsl_vector_set(w->knots, idx++, x_i);

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_col_interp()
  Construct banded ncontrol-by-ncontrol B-spline colocation matrix for interpolation,

X_{ij} = B_j(tau_i)

Inputs: tau - nondecreasing vector, length ncontrol
        XB  - (output) collocation matrix (banded), ncontrol-by-(3*(k-1) + 1) (see below)
        w   - workspace

Return: success/error

Notes:
1) The matrix X is ncontrol-by-ncontrol and banded (p,q) with p = q = k - 1; so the dimensions of the band
storage is ncontrol-by-(3*p + 1)
*/

int
gsl_bspline_col_interp(const gsl_vector * tau,
                       gsl_matrix * XB,
                       gsl_bspline_workspace * w)
{
  const size_t N = w->ncontrol;
  const size_t p = w->spline_order - 1;

  if (tau->size != N)
    {
      GSL_ERROR("tau vector must be length ncontrol", GSL_EBADLEN);
    }
  else if (XB->size1 != N)
    {
      GSL_ERROR("number of rows in XB must match ncontrol", GSL_EBADLEN);
    }
  else if (XB->size2 != 3*p + 1)
    {
      GSL_ERROR("number of columns in XB must be 3*(k-1) + 1", GSL_EBADLEN);
    }
  else
    {
      int status;
      size_t i;

      gsl_matrix_set_zero(XB);

      for (i = 0; i < N; ++i)
        {
          double taui = gsl_vector_get(tau, i);
          gsl_vector_view v;
          size_t istart, idx, j;
          int flag;
          double * ptr;

          idx = gsl_bspline_find_interval (taui, &flag, w);
          if (flag != 0)
            {
              GSL_ERROR("elements of tau must be inside knot interval", GSL_EINVAL);
            }

          j = idx - p; /* column index in full matrix */
          if (i >= j)
            {
              ptr = gsl_matrix_ptr(XB, j, 2*p + i - j);
            }
          else
            {
              GSL_ERROR("runtime error - most likely invalid knot vector", GSL_FAILURE);
            }

          v = gsl_vector_view_array_with_stride(ptr, XB->tda - 1, w->spline_order);

          status = gsl_bspline_basis (taui, &v.vector, &istart, w);
          if (status)
            return status;
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_interp_chermite()
  Compute spline coefficients for Hermite interpolation

Inputs: x  - interpolation sites, length n
        y  - interpolation values y(x), length n
        dy - interpolation values dy/dx(x), length n
        c  - (output) spline coefficients, length 2*(n+2) - k
        w  - workspace

Notes:
1) gsl_bspline_init_hermite() must be called to initialize the knot vector

2) Based on the algorithm in section 4 of:

M. S. Mummy, Hermite interpolation with B-splines, Comp. Aid. Geom. Design, 6, 1989.
*/

int
gsl_bspline_interp_chermite(const gsl_vector * x, const gsl_vector * y,
                            const gsl_vector * dy, gsl_vector * c,
                            const gsl_bspline_workspace * w)
{
  const size_t n = x->size;

  if (n != y->size)
    {
      GSL_ERROR ("x vector does not match y", GSL_EBADLEN);
    }
  else if (n != dy->size)
    {
      GSL_ERROR ("x vector does not match dy", GSL_EBADLEN);
    }
  else if (c->size != 2*n)
    {
      GSL_ERROR ("coefficient vector has wrong size", GSL_EBADLEN);
    }
  else if (w->spline_order != 4) /* this function is specifically for cubic splines */
    {
      GSL_ERROR ("spline_order must be 4", GSL_EDOM);
    }
  else
    {
      double xprev = gsl_vector_get(x, 0);
      size_t i;

      for (i = 0; i < n; ++i)
        {
          double x_i = gsl_vector_get(x, i);
          double y_i = gsl_vector_get(y, i);
          double dy_i = gsl_vector_get(dy, i);
          double * c0 = gsl_vector_ptr(c, 2 * i);
          double * c1 = gsl_vector_ptr(c, 2 * i + 1);

          *c0 = y_i - (x_i - xprev) / 3.0 * dy_i;

          *c1 = y_i;
          if (i < n - 1)
            {
              double xnext = gsl_vector_get(x, i + 1);
              *c1 += (xnext - x_i) / 3.0 * dy_i;
            }

          xprev = x_i;
        }

      return GSL_SUCCESS;
    }
}
