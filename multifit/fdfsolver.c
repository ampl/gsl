/* multifit/fdfsolver.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>

gsl_multifit_fdfsolver *
gsl_multifit_fdfsolver_alloc (const gsl_multifit_fdfsolver_type * T, 
                              size_t n, size_t p)
{
  int status;

  gsl_multifit_fdfsolver * s;

  if (n < p)
    {
      GSL_ERROR_VAL ("insufficient data points, n < p", GSL_EINVAL, 0);
    }

  s = (gsl_multifit_fdfsolver *) calloc (1, sizeof (gsl_multifit_fdfsolver));
  if (s == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for multifit solver struct",
                     GSL_ENOMEM, 0);
    }

  s->x = gsl_vector_calloc (p);

  if (s->x == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
    }

  s->f = gsl_vector_calloc (n);

  if (s->f == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for f", GSL_ENOMEM, 0);
    }

  s->dx = gsl_vector_calloc (p);

  if (s->dx == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for dx", GSL_ENOMEM, 0);
    }

  s->g = gsl_vector_alloc (p);

  if (s->g == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for g", GSL_ENOMEM, 0);
    }

  s->sqrt_wts = gsl_vector_calloc (n);

  if (s->sqrt_wts == 0) 
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for sqrt_wts", GSL_ENOMEM, 0);
    }

  s->state = calloc (1, T->size);

  if (s->state == 0)
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to allocate space for multifit solver state",
                     GSL_ENOMEM, 0);
    }

  s->type = T ;

  status = (s->type->alloc)(s->state, n, p);

  if (status != GSL_SUCCESS)
    {
      gsl_multifit_fdfsolver_free (s);
      GSL_ERROR_VAL ("failed to set solver", status, 0);
    }

  s->fdf = NULL;
  
  s->niter = 0;

  return s;
}

int
gsl_multifit_fdfsolver_set (gsl_multifit_fdfsolver * s, 
                            gsl_multifit_function_fdf * f, 
                            const gsl_vector * x)
{
  return gsl_multifit_fdfsolver_wset(s, f, x, NULL);
}

int
gsl_multifit_fdfsolver_wset (gsl_multifit_fdfsolver * s, 
                             gsl_multifit_function_fdf * f, 
                             const gsl_vector * x,
                             const gsl_vector * wts)
{
  const size_t n = s->f->size;

  if (n != f->n)
    {
      GSL_ERROR ("function size does not match solver", GSL_EBADLEN);
    }
  else if (s->x->size != x->size)
    {
      GSL_ERROR ("vector length does not match solver", GSL_EBADLEN);
    }
  else if (wts != NULL && n != wts->size)
    {
      GSL_ERROR ("weight vector length does not match solver", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      s->fdf = f;
      gsl_vector_memcpy(s->x, x);
      s->niter = 0;

      if (wts)
        {
          for (i = 0; i < n; ++i)
            {
              double wi = gsl_vector_get(wts, i);
              gsl_vector_set(s->sqrt_wts, i, sqrt(wi));
            }
        }
      else
        gsl_vector_set_all(s->sqrt_wts, 1.0);
  
      return (s->type->set) (s->state, s->sqrt_wts, s->fdf, s->x, s->f, s->dx);
    }
}

int
gsl_multifit_fdfsolver_iterate (gsl_multifit_fdfsolver * s)
{
  int status =
    (s->type->iterate) (s->state, s->sqrt_wts, s->fdf, s->x, s->f, s->dx);

  s->niter++;

  return status;
}

/*
gsl_multifit_fdfsolver_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: s - fdfsolver
        maxiter - maximum iterations to allow
        xtol    - tolerance in step x
        gtol    - tolerance in gradient
        ftol    - tolerance in ||f||
        info    - (output) info flag on why iteration terminated
                  1 = stopped due to small step size ||dx|
                  2 = stopped due to small gradient
                  3 = stopped due to small change in f
                  GSL_ETOLX = ||dx|| has converged to within machine
                              precision (and xtol is too small)
                  GSL_ETOLG = ||g||_inf is smaller than machine
                              precision (gtol is too small)
                  GSL_ETOLF = change in ||f|| is smaller than machine
                              precision (ftol is too small)

Return: GSL_SUCCESS if converged, GSL_MAXITER if maxiter exceeded without
converging
*/
int
gsl_multifit_fdfsolver_driver (gsl_multifit_fdfsolver * s,
                               const size_t maxiter,
                               const double xtol,
                               const double gtol,
                               const double ftol,
                               int *info)
{
  int status;
  size_t iter = 0;

  do
    {
      status = gsl_multifit_fdfsolver_iterate (s);

      /*
       * if status is GSL_ENOPROG or GSL_SUCCESS, continue iterating,
       * otherwise the method has converged with a GSL_ETOLx flag
       */
      if (status != GSL_SUCCESS && status != GSL_ENOPROG)
        break;

      /* test for convergence */
      status = gsl_multifit_fdfsolver_test(s, xtol, gtol, ftol, info);
    }
  while (status == GSL_CONTINUE && ++iter < maxiter);

  /*
   * the following error codes mean that the solution has converged
   * to within machine precision, so record the error code in info
   * and return success
   */
  if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
    {
      *info = status;
      status = GSL_SUCCESS;
    }

  /* check if max iterations reached */
  if (iter >= maxiter && status != GSL_SUCCESS)
    status = GSL_EMAXITER;

  return status;
} /* gsl_multifit_fdfsolver_driver() */

int
gsl_multifit_fdfsolver_jac (gsl_multifit_fdfsolver * s, gsl_matrix * J)
{
  const size_t n = s->f->size;
  const size_t p = s->x->size;

  if (n != J->size1 || p != J->size2)
    {
      GSL_ERROR ("Jacobian dimensions do not match solver", GSL_EBADLEN);
    }
  else
    {
      return (s->type->jac) (s->state, J);
    }
} /* gsl_multifit_fdfsolver_jac() */

void
gsl_multifit_fdfsolver_free (gsl_multifit_fdfsolver * s)
{
  RETURN_IF_NULL (s);

  if (s->state)
    {
      (s->type->free) (s->state);
      free (s->state);
    }

  if (s->dx)
    gsl_vector_free (s->dx);

  if (s->x)
    gsl_vector_free (s->x);

  if (s->f)
    gsl_vector_free (s->f);

  if (s->sqrt_wts)
    gsl_vector_free (s->sqrt_wts);

  if (s->g)
    gsl_vector_free (s->g);

  free (s);
}

const char *
gsl_multifit_fdfsolver_name (const gsl_multifit_fdfsolver * s)
{
  return s->type->name;
}

gsl_vector *
gsl_multifit_fdfsolver_position (const gsl_multifit_fdfsolver * s)
{
  return s->x;
}

gsl_vector *
gsl_multifit_fdfsolver_residual (const gsl_multifit_fdfsolver * s)
{
  return s->f;
}

size_t
gsl_multifit_fdfsolver_niter (const gsl_multifit_fdfsolver * s)
{
  return s->niter;
}

/*
gsl_multifit_eval_wf()
  Compute residual vector y with user callback function, and apply
weighting transform if given:

y~ = sqrt(W) y

Inputs: fdf  - callback function
        x    - model parameters
        swts - weight matrix sqrt(W) = sqrt(diag(w1,w2,...,wn))
               set to NULL for unweighted fit
        y    - (output) (weighted) residual vector
               y_i = sqrt(w_i) f_i where f_i is unweighted residual
*/

int
gsl_multifit_eval_wf(gsl_multifit_function_fdf *fdf, const gsl_vector *x,
                     const gsl_vector *swts, gsl_vector *y)
{
  int s = ((*((fdf)->f)) (x, fdf->params, y));
  ++(fdf->nevalf);

  /* y <- sqrt(W) y */
  if (swts)
    gsl_vector_mul(y, swts);

  return s;
}

/*
gsl_multifit_eval_wdf()
  Compute Jacobian matrix J with user callback function, and apply
weighting transform if given:

J~ = sqrt(W) J

Inputs: fdf  - callback function
        x    - model parameters
        swts - weight matrix W = diag(w1,w2,...,wn)
               set to NULL for unweighted fit
        dy   - (output) (weighted) Jacobian matrix
               dy = sqrt(W) dy where dy is unweighted Jacobian
*/

int
gsl_multifit_eval_wdf(gsl_multifit_function_fdf *fdf, const gsl_vector *x,
                      const gsl_vector *swts, gsl_matrix *dy)
{
  int status = ((*((fdf)->df)) (x, fdf->params, dy));

  ++(fdf->nevaldf);

  /* J <- sqrt(W) J */
  if (swts)
    {
      const size_t n = swts->size;
      size_t i;

      for (i = 0; i < n; ++i)
        {
          double swi = gsl_vector_get(swts, i);
          gsl_vector_view v = gsl_matrix_row(dy, i);

          gsl_vector_scale(&v.vector, swi);
        }
    }

  return status;
}
