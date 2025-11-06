/* multifit/gcv.c
 * 
 * Copyright (C) 2016 Patrick Alken
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

/*
 * References:
 *
 * [1] P. C. Hansen, "Discrete Inverse Problems: Insight and Algorithms,"
 * SIAM Press, 2010.
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_min.h>

typedef struct
{
  const gsl_vector * S;
  const gsl_vector * UTy;
  double delta0;
  size_t np;
  gsl_vector * workp;
} gcv_params;

static double gcv_func(double lambda, void * params);

/*
gsl_multifit_linear_gcv_init()
  Initialize Generalized Cross Validation parameters

Inputs: y         - right hand side vector
        reg_param - (output) regularization parameters
        UTy       - (output) U^T y
        delta0    - (output) delta0
        work      - workspace
*/

int
gsl_multifit_linear_gcv_init(const gsl_vector * y,
                             gsl_vector * reg_param,
                             gsl_vector * UTy,
                             double * delta0,
                             gsl_multifit_linear_workspace * work)
{
  const size_t n = y->size;

  if (n != work->n)
    {
      GSL_ERROR("y vector does not match workspace", GSL_EBADLEN);
    }
  else if (UTy->size != work->p)
    {
      GSL_ERROR ("UTy vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      const size_t p = work->p;

      gsl_matrix_view U = gsl_matrix_submatrix(work->A, 0, 0, n, p);
      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);

      const double smax = gsl_vector_get(&S.vector, 0);
      const double smin = gsl_vector_get(&S.vector, p - 1);

      double dr; /* residual error from projection */

      double normy = gsl_blas_dnrm2(y);
      double normUTy;

      /* compute projection UTy = U^T y */
      gsl_blas_dgemv (CblasTrans, 1.0, &U.matrix, y, 0.0, UTy);
      normUTy = gsl_blas_dnrm2(UTy);

      /* dr = ||y||^2 - ||U^T y||^2 */
      dr = (normy + normUTy) * (normy - normUTy);

      /* calculate regularization parameters */
      gsl_multifit_linear_lreg(smin, smax, reg_param);

      if (n > p && dr > 0.0)
        *delta0 = dr;
      else
        *delta0 = 0.0;

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_gcv_curve()
  Calculate Generalized Cross Validation curve for a set
of regularization parameters

Inputs: reg_param - regularization parameters
        UTy       - U^T y vector, size p
        delta0    - delta0
        G         - (output) GCV curve values
        work      - workspace
*/

int
gsl_multifit_linear_gcv_curve(const gsl_vector * reg_param,
                              const gsl_vector * UTy,
                              const double delta0,
                              gsl_vector * G,
                              gsl_multifit_linear_workspace * work)
{
  const size_t n = work->n;
  const size_t p = work->p;
  const size_t N = reg_param->size; /* number of points on GCV curve */

  if (UTy->size != p)
    {
      GSL_ERROR("UTy vector does not match workspace", GSL_EBADLEN);
    }
  else if (G->size != N)
    {
      GSL_ERROR ("size of reg_param and G vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      size_t i;

      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);
      gsl_vector_view workp = gsl_matrix_subcolumn(work->QSI, 0, 0, p);

      gcv_params params;

      params.S = &S.vector;
      params.UTy = UTy;
      params.delta0 = delta0;
      params.np = n - p;
      params.workp = &workp.vector;

      for (i = 0; i < N; ++i)
        {
          double lambdai = gsl_vector_get(reg_param, i);
          double Gi = gcv_func(lambdai, &params);

          gsl_vector_set(G, i, Gi);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_multifit_linear_gcv_min()
  Find regularization parameter which minimizes GCV curve

Inputs: reg_param - regularization parameters
        UTy       - U^T y vector, size p
        G         - GCV curve values
        delta0    - delta0
        lambda    - (output) optimal regularization parameter
        work      - workspace
*/

int
gsl_multifit_linear_gcv_min(const gsl_vector * reg_param,
                            const gsl_vector * UTy,
                            const gsl_vector * G,
                            const double delta0,
                            double * lambda,
                            gsl_multifit_linear_workspace * work)
{
  const size_t n = work->n;
  const size_t p = work->p;
  const size_t npts = reg_param->size; /* number of points on GCV curve */

  if (UTy->size != p)
    {
      GSL_ERROR("UTy vector does not match workspace", GSL_EBADLEN);
    }
  else if (G->size != npts)
    {
      GSL_ERROR ("size of reg_param and G vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int status;
      const size_t max_iter = 500;
      const double tol = 1.0e-4;
      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);
      gsl_vector_view workp = gsl_matrix_subcolumn(work->QSI, 0, 0, p);
      gcv_params params;
      int idxG = (int) gsl_vector_min_index(G);
      double a = gsl_vector_get(reg_param, GSL_MAX(idxG - 1, 0));
      double b = gsl_vector_get(reg_param, GSL_MIN(idxG + 1, (int) npts - 1));
      double m = gsl_vector_get(reg_param, idxG);
      size_t iter = 0;
      gsl_function F;

      /* XXX FIXME */
      gsl_min_fminimizer *min_workspace_p;

      if (idxG == 0 || idxG == ((int)npts - 1))
        {
          /* the minimum is an endpoint of the curve, no need to search */
          *lambda = m;
          return GSL_SUCCESS;
        }

      /* XXX FIXME */
      min_workspace_p = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

      params.S = &S.vector;
      params.UTy = UTy;
      params.delta0 = delta0;
      params.np = n - p;
      params.workp = &workp.vector;

      F.function = gcv_func;
      F.params = &params;

      gsl_min_fminimizer_set(min_workspace_p, &F, m, a, b);

      do
        {
          iter++;
          status = gsl_min_fminimizer_iterate(min_workspace_p);

          a = gsl_min_fminimizer_x_lower(min_workspace_p);
          b = gsl_min_fminimizer_x_upper(min_workspace_p);

          status = gsl_min_test_interval(a, b, 0.0, tol);
        }
      while (status == GSL_CONTINUE && iter < max_iter);

      if (status == GSL_SUCCESS)
        *lambda = gsl_min_fminimizer_minimum(min_workspace_p);
      else
        status = GSL_EMAXITER;

      gsl_min_fminimizer_free(min_workspace_p);

      return status;
    }
}

/*
gsl_multifit_linear_gcv_calc()
  Calculate GCV function G(lambda) for given lambda

Inputs: reg_param - regularization parameters
        UTy       - U^T y vector, size p
        delta0    - delta0
        G         - (output) GCV curve values
        work      - workspace
*/

double
gsl_multifit_linear_gcv_calc(const double lambda,
                             const gsl_vector * UTy,
                             const double delta0,
                             gsl_multifit_linear_workspace * work)
{
  const size_t n = work->n;
  const size_t p = work->p;

  if (UTy->size != p)
    {
      GSL_ERROR_VAL("UTy vector does not match workspace", GSL_EBADLEN, 0.0);
    }
  else
    {
      gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);
      gsl_vector_view workp = gsl_matrix_subcolumn(work->QSI, 0, 0, p);
      gcv_params params;
      double G;

      params.S = &S.vector;
      params.UTy = UTy;
      params.delta0 = delta0;
      params.np = n - p;
      params.workp = &workp.vector;

      G = gcv_func(lambda, &params);

      return G;
    }
}

/*
gsl_multifit_linear_gcv()
  Calculate Generalized Cross Validation curve for a set
of regularization parameters

Inputs: y         - right hand side vector
        reg_param - (output) regularization parameters
        G         - (output) GCV curve values
        lambda    - (output) optimal regularization parameter which
                    minimizes GCV curve
        G_lambda  - (output) G(lambda) value at optimal parameter
        work      - workspace
*/

int
gsl_multifit_linear_gcv(const gsl_vector * y,
                        gsl_vector * reg_param,
                        gsl_vector * G,
                        double * lambda,
                        double * G_lambda,
                        gsl_multifit_linear_workspace * work)
{
  const size_t n = y->size;
  const size_t N = G->size; /* number of points on GCV curve */

  if (n != work->n)
    {
      GSL_ERROR("y vector does not match workspace", GSL_EBADLEN);
    }
  else if (reg_param->size != N)
    {
      GSL_ERROR ("size of reg_param and G vectors do not match",
                 GSL_EBADLEN);
    }
  else
    {
      int status;
      const size_t p = work->p;
      gsl_vector_view UTy = gsl_vector_subvector(work->xt, 0, p);
      double delta0;

      status = gsl_multifit_linear_gcv_init(y, reg_param, &UTy.vector, &delta0, work);
      if (status)
        return status;

      status = gsl_multifit_linear_gcv_curve(reg_param, &UTy.vector, delta0, G, work);
      if (status)
        return status;

      status = gsl_multifit_linear_gcv_min(reg_param, &UTy.vector, G, delta0, lambda, work);
      if (status)
        return status;

      *G_lambda = gsl_multifit_linear_gcv_calc(*lambda, &UTy.vector, delta0, work);

      return GSL_SUCCESS;
    }
}

static double
gcv_func(double lambda, void * params)
{
  gcv_params * par = (gcv_params *) params;
  const gsl_vector *S = par->S;
  const gsl_vector *UTy = par->UTy;
  double delta0 = par->delta0;
  size_t np = par->np;
  gsl_vector *workp = par->workp;
  const size_t p = S->size;
  size_t i;
  double lambda_sq = lambda * lambda;
  double G, d, norm;
  double sumf = 0.0;

  /* compute workp = 1 - filter_factors */
  for (i = 0; i < p; ++i)
    {
      double si = gsl_vector_get(S, i);
      double fi = lambda_sq / (si * si + lambda_sq);
      gsl_vector_set(workp, i, fi);
      sumf += fi;
    }

  d = (double)np + sumf;

  gsl_vector_mul(workp, UTy);
  norm = gsl_blas_dnrm2(workp);

  G = (norm*norm + delta0) / (d * d);

  return G;
}
