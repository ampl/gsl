/* bspline/test_integ.c
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

#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>

static double
poly_func(const double x, void * params)
{
  gsl_vector * coef = (gsl_vector *) params;
  return gsl_poly_eval(coef->data, coef->size, x);
}

int
test_integ_eps(const gsl_function * F, const double a, const double b,
               const double exact, const double tol)
{
  int status = 0;
  gsl_vector * coef = (gsl_vector *) F->params;
  const size_t order = coef->size;
  const size_t ncontrol = 10;
  gsl_bspline_workspace * w = gsl_bspline_alloc_ncontrol(order, ncontrol);
  gsl_vector * y = gsl_vector_alloc(ncontrol);
  gsl_vector * c = gsl_vector_alloc(ncontrol);
  gsl_matrix * G = gsl_matrix_alloc(ncontrol, order);
  double result;

  gsl_bspline_init_uniform(-1.0, 1.0, w);

  gsl_bspline_gram(0, G, w);
  gsl_bspline_proj_rhs(F, y, w);

  gsl_linalg_cholesky_band_decomp(G);
  gsl_linalg_cholesky_band_solve(G, y, c);

  gsl_bspline_calc_integ(a, b, c, &result, w);

  gsl_test_rel(result, exact, tol, "integrate order=%zu a=%g b=%g",
               order, a, b);

  gsl_bspline_free(w);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(G);

  return status;
}

/*
test_integ()
  Test B-spline integration by generating polynomials,
projecting them onto B-spline basis, integrating, and comparing
with Gauss-Legendre result
*/

int
test_integ(gsl_rng * rng_p)
{
  int status = 0;
  const double tol = 1.0e-10;
  size_t order;

  for (order = 1; order < 10; ++order)
    {
      gsl_vector * coef = gsl_vector_alloc(order);
      size_t ngl = order / 2 + 1;
      gsl_integration_glfixed_table * table = gsl_integration_glfixed_table_alloc(ngl);
      double a = 2.0 * gsl_rng_uniform(rng_p) - 1.0; /* integration limits in [-1,1] */
      double b = 2.0 * gsl_rng_uniform(rng_p) - 1.0;
      double lower = GSL_MIN(a, b);
      double upper = GSL_MAX(a, b);
      double result;
      gsl_function F;

      random_vector(-5.0, 5.0, coef, rng_p);

      F.function = poly_func;
      F.params = coef;

      /* compute exact result */
      result = gsl_integration_glfixed(&F, lower, upper, table);

      status += test_integ_eps(&F, lower, upper,  result, tol);
      status += test_integ_eps(&F, upper, lower, -result, tol);

      gsl_vector_free(coef);
      gsl_integration_glfixed_table_free(table);
    }

  return status;
}
