/* multifit_nlinear/fdjac.c
 * 
 * Copyright (C) 2013, 2016 Patrick Alken
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
 * This module contains routines for approximating the Jacobian with
 * finite differences for nonlinear least-squares fitting.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

static int forward_jac(const double h, const gsl_vector *x,
                       const gsl_vector *wts,
                       gsl_multifit_nlinear_fdf *fdf,
                       const gsl_vector *f, gsl_matrix *J);
static int center_jac(const double h, const gsl_vector *x, const gsl_vector *wts,
                      gsl_multifit_nlinear_fdf *fdf, gsl_matrix *J, gsl_vector *work);

/*
forward_jac()
  Compute approximate Jacobian using forward differences

Inputs: h   - finite difference step size
        x   - parameter vector
        wts - data weights
        fdf - fdf struct
        f   - (input) vector of function values f_i(x)
        J   - (output) Jacobian matrix

Return: success or error
*/

static int
forward_jac(const double h, const gsl_vector *x, const gsl_vector *wts,
            gsl_multifit_nlinear_fdf *fdf, const gsl_vector *f, gsl_matrix *J)
{
  int status = 0;
  size_t i, j;
  double delta;

  for (j = 0; j < fdf->p; ++j)
    {
      double xj = gsl_vector_get(x, j);

      /* use column j of J as temporary storage for f(x + dx) */
      gsl_vector_view v = gsl_matrix_column(J, j);

      delta = h * fabs(xj);
      if (delta == 0.0)
        delta = h;

      /* perturb x_j to compute forward difference */
      gsl_vector_set((gsl_vector *) x, j, xj + delta);

      status += gsl_multifit_nlinear_eval_f (fdf, x, wts, &v.vector);
      if (status)
        return status;

      /* restore x_j */
      gsl_vector_set((gsl_vector *) x, j, xj);

      delta = 1.0 / delta;
      for (i = 0; i < fdf->n; ++i)
        {
          double fnext = gsl_vector_get(&v.vector, i);
          double fi = gsl_vector_get(f, i);

          gsl_matrix_set(J, i, j, (fnext - fi) * delta);
        }
    }

  return status;
}

/*
center_jac()
  Compute approximate Jacobian using centered differences

Inputs: h    - finite difference step size
        x    - parameter vector
        wts  - data weights
        fdf  - fdf struct
        J    - (output) Jacobian matrix
        work - additional workspace, size n

Return: success or error
*/

static int
center_jac(const double h, const gsl_vector *x, const gsl_vector *wts,
           gsl_multifit_nlinear_fdf *fdf, gsl_matrix *J, gsl_vector *work)
{
  int status = 0;
  size_t i, j;
  double delta;

  for (j = 0; j < fdf->p; ++j)
    {
      double xj = gsl_vector_get(x, j);

      /* use column j of J as temporary storage for f(x + dx) */
      gsl_vector_view v = gsl_matrix_column(J, j);

      delta = h * fabs(xj);
      if (delta == 0.0)
        delta = h;

      /* perturb x_j to compute forward difference, f(x + 1/2 delta e_j) */
      gsl_vector_set((gsl_vector *) x, j, xj + 0.5 * delta);

      status += gsl_multifit_nlinear_eval_f (fdf, x, wts, &v.vector);
      if (status)
        return status;

      /* perturb x_j to compute backward difference, f(x - 1/2 delta e_j) */
      gsl_vector_set((gsl_vector *) x, j, xj - 0.5 * delta);

      status += gsl_multifit_nlinear_eval_f (fdf, x, wts, work);
      if (status)
        return status;

      /* restore x_j */
      gsl_vector_set((gsl_vector *) x, j, xj);

      delta = 1.0 / delta;
      for (i = 0; i < fdf->n; ++i)
        {
          double fnext = gsl_vector_get(&v.vector, i);
          double fprev = gsl_vector_get(work, i);

          gsl_matrix_set(J, i, j, (fnext - fprev) * delta);
        }
    }

  return status;
}

/*
gsl_multifit_nlinear_df()
  Compute approximate Jacobian using finite differences

Inputs: h      - finite difference step size
        fdtype - finite difference method
        x      - parameter vector
        wts    - data weights (set to NULL if not needed)
        fdf    - fdf
        f      - (input) function values f_i(x)
        J      - (output) approximate (weighted) Jacobian matrix, sqrt(W) * J
        work   - additional workspace for centered differences, size n

Return: success or error
*/

int
gsl_multifit_nlinear_df(const double h, const gsl_multifit_nlinear_fdtype fdtype,
                        const gsl_vector *x, const gsl_vector *wts,
                        gsl_multifit_nlinear_fdf *fdf,
                        const gsl_vector *f, gsl_matrix *J, gsl_vector *work)
{
  int status;

  if (fdtype == GSL_MULTIFIT_NLINEAR_FWDIFF)
    {
      status = forward_jac(h, x, wts, fdf, f, J);
    }
  else if (fdtype == GSL_MULTIFIT_NLINEAR_CTRDIFF)
    {
      status = center_jac(h, x, wts, fdf, J, work);
    }
  else
    {
      GSL_ERROR("invalid specified fdtype", GSL_EINVAL);
    }

  return status;
}
