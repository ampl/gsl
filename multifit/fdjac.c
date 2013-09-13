/* multifit/fdjac.c
 * 
 * Copyright (C) 2013 Patrick Alken
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
 *
 *
 * This module contains routines for approximating the Jacobian with finite
 * differences for nonlinear least-squares fitting.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

static int fdjac(const gsl_vector *x, gsl_multifit_function_fdf *fdf,
                 const gsl_vector *f, gsl_matrix *J);

/*
fdjac()
  Compute approximate Jacobian using forward differences

Inputs: x   - parameter vector
        fdf - fdf struct
        f   - (input) vector of function values f_i(x)
        J   - (output) Jacobian matrix

Return: success or error
*/

static int
fdjac(const gsl_vector *x, gsl_multifit_function_fdf *fdf,
      const gsl_vector *f, gsl_matrix *J)
{
  int status = 0;
  size_t i, j;
  double h;
  const double epsfcn = 0.0;
  double eps = sqrt(GSL_MAX(epsfcn, GSL_DBL_EPSILON));

  for (j = 0; j < fdf->p; ++j)
    {
      double xj = gsl_vector_get(x, j);

      /* use column j of J as temporary storage for f(x + dx) */
      gsl_vector_view v = gsl_matrix_column(J, j);

      h = eps * fabs(xj);
      if (h == 0.0)
        h = eps;

      /* perturb x_j to compute forward difference */
      gsl_vector_set((gsl_vector *) x, j, xj + h);

      status += GSL_MULTIFIT_FN_EVAL_F (fdf, x, &v.vector);
      if (status)
        return status;

      /* restore x_j */
      gsl_vector_set((gsl_vector *) x, j, xj);

      h = 1.0 / h;
      for (i = 0; i < fdf->n; ++i)
        {
          double fnext = gsl_vector_get(&v.vector, i);
          double fi = gsl_vector_get(f, i);

          gsl_matrix_set(J, i, j, (fnext - fi) * h);
        }
    }

  return status;
} /* fdjac() */

/*
gsl_multifit_fdfsolver_dif_df()
  Compute approximate Jacobian using finite differences

Inputs: x   - parameter vector
        fdf - fdf
        f   - (input) function values f_i(x)
        J   - (output) approximate Jacobian matrix

Return: success or error
*/

int
gsl_multifit_fdfsolver_dif_df(const gsl_vector *x, gsl_multifit_function_fdf *fdf,
                              const gsl_vector *f, gsl_matrix *J)
{
  return fdjac(x, fdf, f, J);
} /* gsl_multifit_fdfsolver_dif_df() */

/*
gsl_multifit_fdfsolver_dif_fdf()
  Compute function values (analytic) and approximate Jacobian using finite
differences

Inputs: x      - parameter vector
        fdf    - fdf
        f      - (output) function values f_i(x)
        J      - (output) approximate Jacobian matrix

Return: success or error
*/

int
gsl_multifit_fdfsolver_dif_fdf(const gsl_vector *x, gsl_multifit_function_fdf *fdf,
                               gsl_vector *f, gsl_matrix *J)
{
  int status = 0;

  status = GSL_MULTIFIT_FN_EVAL_F(fdf, x, f);
  if (status)
    return status;

  status = fdjac(x, fdf, f, J);
  if (status)
    return status;

  return status;
} /* gsl_multifit_fdfsolver_dif_fdf() */
