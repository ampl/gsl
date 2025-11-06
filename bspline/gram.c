/* bspline/gram.c
 *
 * Copyright (C) 2018, 2019, 2020 Patrick Alken
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

/*
gsl_bspline_gram()
  Compute the Gram matrix of integrals of products of B-splines.

G_{ij} = \int_{min(knots)}^{max(knots)} N^{(nderiv)}_i(x) N^{(nderiv)}_j(x) dx

The matrix is ncontrol-by-ncontrol, symmetric and banded, with lower bandwidth k-1.

Inputs: nderiv - derivative order
        G      - (output) Gram matrix in symmetric banded format,
                 ncontrol-by-order
        w      - workspace

Notes:
1) See Algorithm 5.22 [pg. 203] of L. Schumaker, "Spline Functions: Basic Theory", 3rd ed.
2) The algorithm uses Gauss-Legendre quadrature of size spline_order on each knot interval
*/

int
gsl_bspline_gram(const size_t nderiv, gsl_matrix * G, gsl_bspline_workspace * w)
{
  const double a = gsl_vector_get(w->knots, 0);
  const double b = gsl_vector_get(w->knots, w->knots->size - 1);

  return gsl_bspline_gram_interval(a, b, nderiv, G, w);
}

/*
gsl_bspline_gram_interval()
  Compute the Gram matrix of integrals of products of B-splines on the
interval [a,b]

G_{ij} = \int_a^b N^{(nderiv)}_i(x) N^{(nderiv)}_j(x) dx

The matrix is ncontrol-by-ncontrol, symmetric and banded, with lower bandwidth k-1.

Inputs: a      - lower integration limit
        b      - upper integration limit
        nderiv - derivative order
        G      - (output) Gram matrix in symmetric banded format,
                 ncontrol-by-order
        w      - workspace

Notes:
1) See Algorithm 5.22 [pg. 203] of L. Schumaker, "Spline Functions: Basic Theory", 3rd ed.
2) The algorithm uses Gauss-Legendre quadrature of size spline_order on each knot interval
*/

int
gsl_bspline_gram_interval(const double a, const double b, const size_t nderiv,
                          gsl_matrix * G, gsl_bspline_workspace * w)
{
  const size_t ncontrol = w->ncontrol;
  const size_t order = w->spline_order;

  if (b < a)
    {
      GSL_ERROR ("b must be greater than a", GSL_EDOM);
    }
  else if (G->size1 != ncontrol)
    {
      GSL_ERROR ("first matrix dimension must equal ncontrol", GSL_EBADLEN);
    }
  else if (G->size2 != order)
    {
      GSL_ERROR ("second matrix dimension must equal spline order", GSL_EBADLEN);
    }
  else if (nderiv >= order)
    {
      /* quick return */
      gsl_matrix_set_zero(G);
      return GSL_SUCCESS;
    }
  else
    {
      const size_t ngl = order - nderiv; /* number of Gauss-Legendre points */
      const size_t m = (ngl + 1) >> 1;
      size_t nu;
      gsl_integration_glfixed_table * gltable;
      gsl_matrix_view N1 = gsl_matrix_submatrix(w->dB, 0, 0, order, nderiv + 1);
      gsl_matrix_view N2 = gsl_matrix_submatrix(w->dB, 0, nderiv + 1, order, nderiv + 1);
      size_t ustart = 0;
      size_t istart, i, j;

      gsl_matrix_set_zero(G);

      /* quick return */
      if (a == b)
        return GSL_SUCCESS;

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
              double t0 = GSL_MAX(gsl_vector_get(w->knots, nu), a);
              double t1 = GSL_MIN(gsl_vector_get(w->knots, nu + 1), b);

              if (t1 > t0)
                {
                  double A = 0.5 * (t1 - t0);
                  double B = t0 + A; /* 0.5 * (t1 + t0) */
                  double w0 = A * gltable->w[0];

                  gsl_bspline_basis_deriv(B, nderiv, &N1.matrix, &istart, w);

                  for (i = 0; i < order; ++i)
                    {
                      double Ni = gsl_matrix_get(&N1.matrix, i, nderiv);

                      for (j = 0; j <= i; ++j)
                        {
                          double Nj = gsl_matrix_get(&N1.matrix, j, nderiv);
                          double * Gij = gsl_matrix_ptr(G, istart + j, i - j);
                          *Gij += w0 * (Ni * Nj);
                        }
                    }
                }
            }
        }

      for (nu = 0; nu < ncontrol + order - 1; ++nu)
        {
          double t0 = GSL_MAX(gsl_vector_get(w->knots, nu), a);
          double t1 = GSL_MIN(gsl_vector_get(w->knots, nu + 1), b);

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
                  gsl_bspline_basis_deriv(B + Ax, nderiv, &N1.matrix, &istart, w);
                  gsl_bspline_basis_deriv(B - Ax, nderiv, &N2.matrix, &istart, w);

                  for (i = 0; i < order; ++i)
                    {
                      double N1i = gsl_matrix_get(&N1.matrix, i, nderiv);
                      double N2i = gsl_matrix_get(&N2.matrix, i, nderiv);

                      for (j = 0; j <= i; ++j)
                        {
                          double N1j = gsl_matrix_get(&N1.matrix, j, nderiv);
                          double N2j = gsl_matrix_get(&N2.matrix, j, nderiv);
                          double * Gij = gsl_matrix_ptr(G, istart + j, i - j);
                          *Gij += wu * (N1i * N1j + N2i * N2j);
                        }
                    }
                }
            }
        }

      gsl_integration_glfixed_table_free(gltable);

      return GSL_SUCCESS;
    }
}

/*
gsl_bspline_oprod()
  Compute the outer product matrix of B-spline basis.

A_{ij}(x) = N^{(nderiv)}_i(x) N^{(nderiv)}_j(x)

The matrix is ncontrol-by-ncontrol, symmetric and banded, with lower bandwidth k-1.

Inputs: nderiv - derivative order
        x      - point to evaluate the B-spline basis functions
        A      - (output) outer product matrix in symmetric banded format,
                 ncontrol-by-order
        w      - workspace
*/

int
gsl_bspline_oprod(const size_t nderiv, const double x, gsl_matrix * A, gsl_bspline_workspace * w)
{
  const size_t ncontrol = w->ncontrol;
  const size_t order = w->spline_order;

  if (A->size1 != ncontrol)
    {
      GSL_ERROR ("first matrix dimension must equal ncontrol", GSL_EBADLEN);
    }
  else if (A->size2 != order)
    {
      GSL_ERROR ("second matrix dimension must equal spline order", GSL_EBADLEN);
    }
  else if (nderiv >= order)
    {
      /* quick return */
      gsl_matrix_set_zero(A);
      return GSL_SUCCESS;
    }
  else
    {
      gsl_vector_const_view N = gsl_matrix_const_column(w->dB, nderiv);
      size_t i, j, istart;

      gsl_matrix_set_zero(A);

      gsl_bspline_basis_deriv(x, nderiv, w->dB, &istart, w);

      for (i = 0; i < order; ++i)
        {
          double Ni = gsl_vector_get(&N.vector, i);

          for (j = 0; j <= i; ++j)
            {
              double Nj = gsl_vector_get(&N.vector, j);
              gsl_matrix_set(A, istart + j, i - j, Ni * Nj);
            }
        }

      return GSL_SUCCESS;
    }
}
