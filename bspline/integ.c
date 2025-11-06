/* bspline/integ.c
 *
 * Copyright (C) 2020 Patrick Alken
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

/* routines related to evaluation of B-splines integrals */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

/*
gsl_bspline_calc_integ()
  Compute

result = \int_a^b f(x) dx = \sum_i c_i \int_a^b B_i(x) dx

Inputs: a      - lower integration limit
        b      - upper integration limit
        c      - spline coefficients, length ncontrol
        result - (output) integral value
        w      - workspace

Return: success/error
*/

int
gsl_bspline_calc_integ(const double a, const double b,
                       const gsl_vector * c, double * result,
                       gsl_bspline_workspace * w)
{
  if (c->size != w->ncontrol)
    {
      GSL_ERROR("c vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status;
      gsl_vector_view bint = gsl_vector_subvector(w->work, 0, w->ncontrol);

      status = gsl_bspline_basis_integ(a, b, &bint.vector, w);
      if (status)
        return status;

      gsl_blas_ddot(c, &bint.vector, result);

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_basis_integ()
  Compute integrals of B-spline basis functions over [a,b]

Inputs: a    - lower integration limit
        b    - upper integration limit
        bint - (output) vector of B-spline integrals, length ncontrol
               bint(i) = \int_a^b B_i(x) dx
        w    - workspace

Return: success/error
*/

int
gsl_bspline_basis_integ(const double a, const double b,
                        gsl_vector * bint,
                        gsl_bspline_workspace * w)
{
  if (bint->size != w->ncontrol)
    {
      GSL_ERROR("bint vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      const size_t order = w->spline_order;
      const size_t ncontrol = w->ncontrol;
      const size_t ngl = order / 2 + 1; /* number of Gauss-Legendre points */
      const size_t m = (ngl + 1) >> 1;
      size_t nu;
      gsl_integration_glfixed_table * gltable;
      gsl_vector_view N1 = gsl_matrix_column(w->dB, 0);
      gsl_vector_view N2 = gsl_matrix_column(w->dB, 1);
      double lower, upper;
      size_t ustart = 0;
      size_t istart, i;

      gsl_vector_set_zero(bint);

      /* quick return */
      if (a == b)
        return GSL_SUCCESS;

      lower = GSL_MIN(a, b);
      upper = GSL_MAX(a, b);

      /* for reasonably small spline orders, no dynamic allocation occurs */
      gltable = gsl_integration_glfixed_table_alloc(ngl);
      if (gltable == NULL)
        {
          GSL_ERROR ("failed to initialize Gauss-Legendre quadrature points", GSL_EFAILED);
        }

      if (ngl & 1)
        {
          /* odd number of Gauss-Legendre points; initialize output matrix
           * with results of middle (w0) weight */

          ustart = 1;

          for (nu = 0; nu < ncontrol + order - 1; ++nu)
            {
              double t0 = GSL_MAX(gsl_vector_get(w->knots, nu), lower);
              double t1 = GSL_MIN(gsl_vector_get(w->knots, nu + 1), upper);

              if (t1 > t0)
                {
                  double A = 0.5 * (t1 - t0);
                  double B = t0 + A; /* 0.5 * (t1 + t0) */
                  double w0 = A * gltable->w[0];

                  gsl_bspline_basis(B, &N1.vector, &istart, w);

                  for (i = 0; i < order; ++i)
                    {
                      double Ni = gsl_vector_get(&N1.vector, i);
                      double * ptr = gsl_vector_ptr(bint, istart + i);
                      *ptr += w0 * Ni;
                    }
                }
            }
        }

      for (nu = 0; nu < ncontrol + order - 1; ++nu)
        {
          double t0 = GSL_MAX(gsl_vector_get(w->knots, nu), lower);
          double t1 = GSL_MIN(gsl_vector_get(w->knots, nu + 1), upper);

          if (t1 > t0)
            {
              double A = 0.5 * (t1 - t0);
              double B = t0 + A; /* 0.5 * (t1 + t0) */
              size_t u;

              for (u = ustart; u < m; ++u)
                {
                  double Ax = A * gltable->x[u];
                  double wu = A * gltable->w[u];

                  /* both of these calls will return the same istart, since we are inside
                   * the knot interval [t0, t1] */
                  gsl_bspline_basis(B + Ax, &N1.vector, &istart, w);
                  gsl_bspline_basis(B - Ax, &N2.vector, &istart, w);

                  for (i = 0; i < order; ++i)
                    {
                      double N1i = gsl_vector_get(&N1.vector, i);
                      double N2i = gsl_vector_get(&N2.vector, i);
                      double * ptr = gsl_vector_ptr(bint, istart + i);
                      *ptr += wu * (N1i + N2i);
                    }
                }
            }
        }

      if (b < a)
        gsl_vector_scale(bint, -1.0);

      gsl_integration_glfixed_table_free(gltable);

      return GSL_SUCCESS;
    }
}
