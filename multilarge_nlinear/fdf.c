/* multilarge_nlinear/fdf.c
 * 
 * Copyright (C) 2015, 2016 Patrick Alken
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multilarge_nlinear.h>

gsl_multilarge_nlinear_workspace *
gsl_multilarge_nlinear_alloc (const gsl_multilarge_nlinear_type * T, 
                              const gsl_multilarge_nlinear_parameters * params,
                              const size_t n, const size_t p)
{
  gsl_multilarge_nlinear_workspace * w;

  if (n < p)
    {
      GSL_ERROR_VAL ("insufficient data points, n < p", GSL_EINVAL, 0);
    }

  w = calloc (1, sizeof (gsl_multilarge_nlinear_workspace));
  if (w == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for workspace",
                     GSL_ENOMEM, 0);
    }

  w->n = n;
  w->p = p;
  w->type = T;
  w->fdf = NULL;
  w->niter = 0;
  w->params = *params;

  /* the cgst method uses its own built-in linear solver */
  if (w->params.trs == gsl_multilarge_nlinear_trs_cgst)
    {
      w->params.solver = gsl_multilarge_nlinear_solver_none;
    }

  w->x = gsl_vector_calloc (p);
  if (w->x == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
    }

  w->f = gsl_vector_calloc (n);
  if (w->f == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for f", GSL_ENOMEM, 0);
    }

  w->dx = gsl_vector_calloc (p);
  if (w->dx == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for dx", GSL_ENOMEM, 0);
    }

  w->g = gsl_vector_alloc (p);
  if (w->g == 0) 
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for g", GSL_ENOMEM, 0);
    }

  if (w->params.solver == gsl_multilarge_nlinear_solver_cholesky ||
      w->params.solver == gsl_multilarge_nlinear_solver_mcholesky)
    {
      w->JTJ = gsl_matrix_alloc (p, p);
      if (w->JTJ == 0) 
        {
          gsl_multilarge_nlinear_free (w);
          GSL_ERROR_VAL ("failed to allocate space for JTJ", GSL_ENOMEM, 0);
        }
    }

  w->sqrt_wts_work = gsl_vector_calloc (n);
  if (w->sqrt_wts_work == 0)
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for weights", GSL_ENOMEM, 0);
    }

  w->state = (T->alloc)(&(w->params), n, p);
  if (w->state == 0)
    {
      gsl_multilarge_nlinear_free (w);
      GSL_ERROR_VAL ("failed to allocate space for state", GSL_ENOMEM, 0);
    }

  return w;
}

void
gsl_multilarge_nlinear_free (gsl_multilarge_nlinear_workspace * w)
{
  RETURN_IF_NULL (w);

  if (w->state)
    (w->type->free) (w->state);

  if (w->dx)
    gsl_vector_free (w->dx);

  if (w->x)
    gsl_vector_free (w->x);

  if (w->f)
    gsl_vector_free (w->f);

  if (w->sqrt_wts_work)
    gsl_vector_free (w->sqrt_wts_work);

  if (w->g)
    gsl_vector_free (w->g);

  if (w->JTJ)
    gsl_matrix_free (w->JTJ);

  free (w);
}

gsl_multilarge_nlinear_parameters
gsl_multilarge_nlinear_default_parameters(void)
{
  gsl_multilarge_nlinear_parameters params;

  params.trs = gsl_multilarge_nlinear_trs_lm;
  params.scale = gsl_multilarge_nlinear_scale_more;
  params.solver = gsl_multilarge_nlinear_solver_cholesky;
  params.fdtype = GSL_MULTILARGE_NLINEAR_FWDIFF;
  params.factor_up = 3.0;
  params.factor_down = 2.0;
  params.avmax = 0.75;
  params.h_df = GSL_SQRT_DBL_EPSILON;
  params.h_fvv = 0.01;
  params.max_iter = 0;
  params.tol = 1.0e-6;

  return params;
}

int
gsl_multilarge_nlinear_init (const gsl_vector * x,
                             gsl_multilarge_nlinear_fdf * fdf,
                             gsl_multilarge_nlinear_workspace * w)
{
  return gsl_multilarge_nlinear_winit(x, NULL, fdf, w);
}

int
gsl_multilarge_nlinear_winit (const gsl_vector * x,
                              const gsl_vector * wts,
                              gsl_multilarge_nlinear_fdf * fdf, 
                              gsl_multilarge_nlinear_workspace * w)
{
  const size_t n = w->f->size;

  if (n != fdf->n)
    {
      GSL_ERROR ("function size does not match workspace", GSL_EBADLEN);
    }
  else if (w->x->size != x->size)
    {
      GSL_ERROR ("vector length does not match workspace", GSL_EBADLEN);
    }
  else if (wts != NULL && n != wts->size)
    {
      GSL_ERROR ("weight vector length does not match workspace", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* initialize counters for function and Jacobian evaluations */
      fdf->nevalf = 0;
      fdf->nevaldfu = 0;
      fdf->nevaldf2 = 0;
      fdf->nevalfvv = 0;

      w->fdf = fdf;
      gsl_vector_memcpy(w->x, x);
      w->niter = 0;

      if (wts)
        {
          w->sqrt_wts = w->sqrt_wts_work;

          for (i = 0; i < n; ++i)
            {
              double wi = gsl_vector_get(wts, i);
              gsl_vector_set(w->sqrt_wts, i, sqrt(wi));
            }
        }
      else
        {
          w->sqrt_wts = NULL;
        }
  
      return (w->type->init) (w->state, w->sqrt_wts, w->fdf,
                              w->x, w->f, w->g, w->JTJ);
    }
}

int
gsl_multilarge_nlinear_iterate (gsl_multilarge_nlinear_workspace * w)
{
  int status =
    (w->type->iterate) (w->state, w->sqrt_wts, w->fdf,
                        w->x, w->f, w->g, w->JTJ, w->dx);

  w->niter++;

  return status;
}

double
gsl_multilarge_nlinear_avratio (const gsl_multilarge_nlinear_workspace * w)
{
  return (w->type->avratio) (w->state);
}

int
gsl_multilarge_nlinear_rcond (double * rcond, const gsl_multilarge_nlinear_workspace * w)
{
  int status = (w->type->rcond) (rcond, w->JTJ, w->state);
  return status;
}

int
gsl_multilarge_nlinear_covar (gsl_matrix * covar, gsl_multilarge_nlinear_workspace * w)
{
  if (covar->size1 != covar->size2)
    {
      GSL_ERROR ("covariance matrix must be square", GSL_ENOTSQR);
    }
  else if (covar->size1 != w->p)
    {
      GSL_ERROR ("covariance matrix does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status = (w->type->covar) (w->JTJ, covar, w->state);
      return status;
    }
}

/*
gsl_multilarge_nlinear_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: maxiter  - maximum iterations to allow
        xtol     - tolerance in step x
        gtol     - tolerance in gradient
        ftol     - tolerance in ||f||
        callback - callback function to call each iteration
        callback_params - parameters to pass to callback function
        info     - (output) info flag on why iteration terminated
                   1 = stopped due to small step size ||dx|
                   2 = stopped due to small gradient
                   3 = stopped due to small change in f
                   GSL_ETOLX = ||dx|| has converged to within machine
                               precision (and xtol is too small)
                   GSL_ETOLG = ||g||_inf is smaller than machine
                               precision (gtol is too small)
                   GSL_ETOLF = change in ||f|| is smaller than machine
                               precision (ftol is too small)
        w        - workspace

Return:
GSL_SUCCESS if converged
GSL_MAXITER if maxiter exceeded without converging
GSL_ENOPROG if no accepted step found on first iteration
*/

int
gsl_multilarge_nlinear_driver (const size_t maxiter,
                               const double xtol,
                               const double gtol,
                               const double ftol,
                               void (*callback)(const size_t iter, void *params,
                                                const gsl_multilarge_nlinear_workspace *w),
                               void *callback_params,
                               int *info,
                               gsl_multilarge_nlinear_workspace * w)
{
  int status;
  size_t iter = 0;

  /* call user callback function prior to any iterations
   * with initial system state */
  if (callback)
    callback(iter, callback_params, w);

  do
    {
      status = gsl_multilarge_nlinear_iterate (w);

      /*
       * If the solver reports no progress on the first iteration,
       * then it didn't find a single step to reduce the
       * cost function and more iterations won't help so return.
       *
       * If we get a no progress flag on subsequent iterations,
       * it means we did find a good step in a previous iteration,
       * so continue iterating since the solver has now reset
       * mu to its initial value.
       */
      if (status == GSL_ENOPROG && iter == 0)
        {
          *info = status;
          return GSL_EMAXITER;
        }

      ++iter;

      if (callback)
        callback(iter, callback_params, w);

      /* test for convergence */
      status = gsl_multilarge_nlinear_test(xtol, gtol, ftol, info, w);
    }
  while (status == GSL_CONTINUE && iter < maxiter);

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
} /* gsl_multilarge_nlinear_driver() */

const char *
gsl_multilarge_nlinear_name (const gsl_multilarge_nlinear_workspace * w)
{
  return w->type->name;
}

gsl_vector *
gsl_multilarge_nlinear_position (const gsl_multilarge_nlinear_workspace * w)
{
  return w->x;
}

gsl_vector *
gsl_multilarge_nlinear_residual (const gsl_multilarge_nlinear_workspace * w)
{
  return w->f;
}

gsl_vector *
gsl_multilarge_nlinear_step (const gsl_multilarge_nlinear_workspace * w)
{
  return w->dx;
}

size_t
gsl_multilarge_nlinear_niter (const gsl_multilarge_nlinear_workspace * w)
{
  return w->niter;
}

const char *
gsl_multilarge_nlinear_trs_name (const gsl_multilarge_nlinear_workspace * w)
{
  return w->params.trs->name;
}

/*
gsl_multilarge_nlinear_eval_f()
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
gsl_multilarge_nlinear_eval_f(gsl_multilarge_nlinear_fdf *fdf,
                              const gsl_vector *x,
                              const gsl_vector *swts,
                              gsl_vector *y)
{
  int s = ((*((fdf)->f)) (x, fdf->params, y));

  ++(fdf->nevalf);

  /* y <- sqrt(W) y */
  if (swts)
    gsl_vector_mul(y, swts);

  return s;
}

/*
gsl_multilarge_nlinear_eval_df()
  Compute Jacobian matrix-vector product:
  
v = J * u

or

v = J^T u

Inputs: TransJ - use J or J^T
        x      - model parameters
        f      - residual vector f(x)
        u      - input vector u
        swts   - weight matrix W = diag(w1,w2,...,wn)
                 set to NULL for unweighted fit
        h      - finite difference step size
        fdtype - finite difference method
        fdf    - callback function
        v      - (output) vector v
        JTJ    - (output) matrix J^T J
        work   - workspace for finite difference, size n
*/

int
gsl_multilarge_nlinear_eval_df(const CBLAS_TRANSPOSE_t TransJ,
                               const gsl_vector *x,
                               const gsl_vector *f,
                               const gsl_vector *u,
                               const gsl_vector *swts,
                               const double h,
                               const gsl_multilarge_nlinear_fdtype fdtype,
                               gsl_multilarge_nlinear_fdf *fdf,
                               gsl_vector *v,
                               gsl_matrix *JTJ,
                               gsl_vector *work)
{
  const size_t n = fdf->n;
  const size_t p = fdf->p;

  if (u != NULL && ((TransJ == CblasNoTrans && u->size != p) ||
                    (TransJ == CblasTrans && u->size != n)))
    {
      GSL_ERROR("u vector has wrong size", GSL_EBADLEN);
    }
  else if (v != NULL && ((TransJ == CblasNoTrans && v->size != n) ||
                         (TransJ == CblasTrans && v->size != p)))
    {
      GSL_ERROR("v vector has wrong size", GSL_EBADLEN);
    }
  else if (JTJ != NULL && ((JTJ->size1 != p) || (JTJ->size2 != p)))
    {
      GSL_ERROR("JTJ matrix has wrong size", GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;

      if (fdf->df)
        {
          /* call user-supplied function */
          status = ((*((fdf)->df)) (TransJ, x, u, fdf->params, v, JTJ));

          if (v)
            ++(fdf->nevaldfu);

          if (JTJ)
            ++(fdf->nevaldf2);
        }
      else
        {
#if 0
          /* use finite difference Jacobian approximation */
          status = gsl_multilarge_nlinear_df(h, fdtype, x, swts, fdf, f, df, work);
#endif
        }

      return status;
    }
}

/*
gsl_multilarge_nlinear_eval_fvv()
  Compute second direction derivative vector yvv with user
callback function, and apply weighting transform if given:

yvv~ = sqrt(W) yvv

Inputs: h    - step size for finite difference, if needed
        x    - model parameters, size p
        v    - unscaled geodesic velocity vector, size p
        f    - residual vector f(x), size n
        swts - weight matrix sqrt(W) = sqrt(diag(w1,w2,...,wn))
               set to NULL for unweighted fit
        fdf  - callback function
        yvv  - (output) (weighted) second directional derivative vector
               yvv_i = sqrt(w_i) fvv_i where f_i is unweighted
        work - workspace, size p
*/

int
gsl_multilarge_nlinear_eval_fvv(const double h,
                                const gsl_vector *x,
                                const gsl_vector *v,
                                const gsl_vector *f,
                                const gsl_vector *swts,
                                gsl_multilarge_nlinear_fdf *fdf,
                                gsl_vector *yvv,
                                gsl_vector *work)
{
  int status = GSL_SUCCESS;
  
  if (fdf->fvv != NULL)
    {
      /* call user-supplied function */
      status = ((*((fdf)->fvv)) (x, v, fdf->params, yvv));
      ++(fdf->nevalfvv);
    }
  else
    {
#if 0
      /* use finite difference approximation */
      status = gsl_multilarge_nlinear_fdfvv(h, x, v, f, J,
                                          swts, fdf, yvv, work);
#endif
    }

  /* yvv <- sqrt(W) yvv */
  if (swts)
    gsl_vector_mul(yvv, swts);

  return status;
}
