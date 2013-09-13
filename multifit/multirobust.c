/* multirobust.c
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
 * This module contains routines related to robust linear least squares. The
 * algorithm used closely follows the publications:
 *
 * [1] DuMouchel, W. and F. O'Brien (1989), "Integrating a robust
 * option into a multiple regression computing environment,"
 * Computer Science and Statistics:  Proceedings of the 21st
 * Symposium on the Interface, American Statistical Association
 *
 * [2] Street, J.O., R.J. Carroll, and D. Ruppert (1988), "A note on
 * computing robust regression estimates via iteratively
 * reweighted least squares," The American Statistician, v. 42, 
 * pp. 152-154.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static int robust_test_convergence(const gsl_vector *c_prev, const gsl_vector *c,
                                   const double tol);
static double robust_madsigma(const gsl_vector *x, gsl_multifit_robust_workspace *w);
static double robust_robsigma(const gsl_vector *r, const double s,
                              const double tune, gsl_multifit_robust_workspace *w);
static double robust_sigma(const double s_ols, const double s_rob,
                           gsl_multifit_robust_workspace *w);
static int robust_covariance(const double sigma, gsl_matrix *cov,
                             gsl_multifit_robust_workspace *w);

/*
gsl_multifit_robust_alloc
  Allocate a robust workspace

Inputs: T - robust weighting algorithm
        n - number of observations
        p - number of model parameters

Return: pointer to workspace
*/

gsl_multifit_robust_workspace *
gsl_multifit_robust_alloc(const gsl_multifit_robust_type *T,
                          const size_t n, const size_t p)
{
  gsl_multifit_robust_workspace *w;

  if (n < p)
    {
      GSL_ERROR_VAL("observations n must be >= p", GSL_EINVAL, 0);
    }

  w = calloc(1, sizeof(gsl_multifit_robust_workspace));
  if (w == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for multifit_robust struct",
                    GSL_ENOMEM, 0);
    }

  w->n = n;
  w->p = p;
  w->type = T;
  w->maxiter = 100; /* maximum iterations */
  w->tune = w->type->tuning_default;

  w->multifit_p = gsl_multifit_linear_alloc(n, p);
  if (w->multifit_p == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for multifit_linear struct",
                    GSL_ENOMEM, 0);
    }

  w->r = gsl_vector_alloc(n);
  if (w->r == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for residuals",
                    GSL_ENOMEM, 0);
    }

  w->weights = gsl_vector_alloc(n);
  if (w->weights == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for weights", GSL_ENOMEM, 0);
    }

  w->c_prev = gsl_vector_alloc(p);
  if (w->c_prev == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for c_prev", GSL_ENOMEM, 0);
    }

  w->resfac = gsl_vector_alloc(n);
  if (w->resfac == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for residual factors",
                    GSL_ENOMEM, 0);
    }

  w->psi = gsl_vector_alloc(n);
  if (w->psi == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for psi", GSL_ENOMEM, 0);
    }

  w->dpsi = gsl_vector_alloc(n);
  if (w->dpsi == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for dpsi", GSL_ENOMEM, 0);
    }

  w->QSI = gsl_matrix_alloc(p, p);
  if (w->QSI == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for QSI", GSL_ENOMEM, 0);
    }

  w->D = gsl_vector_alloc(p);
  if (w->D == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for D", GSL_ENOMEM, 0);
    }

  w->workn = gsl_vector_alloc(n);
  if (w->workn == 0)
    {
      GSL_ERROR_VAL("failed to allocate space for workn", GSL_ENOMEM, 0);
    }

  w->stats.sigma_ols = 0.0;
  w->stats.sigma_mad = 0.0;
  w->stats.sigma_rob = 0.0;
  w->stats.sigma = 0.0;
  w->stats.Rsq = 0.0;
  w->stats.adj_Rsq = 0.0;
  w->stats.rmse = 0.0;
  w->stats.sse = 0.0;
  w->stats.dof = n - p;
  w->stats.weights = w->weights;
  w->stats.r = w->r;

  return w;
} /* gsl_multifit_robust_alloc() */

/*
gsl_multifit_robust_free()
  Free memory associated with robust workspace
*/

void
gsl_multifit_robust_free(gsl_multifit_robust_workspace *w)
{
  if (w->multifit_p)
    gsl_multifit_linear_free(w->multifit_p);

  if (w->r)
    gsl_vector_free(w->r);

  if (w->weights)
    gsl_vector_free(w->weights);

  if (w->c_prev)
    gsl_vector_free(w->c_prev);

  if (w->resfac)
    gsl_vector_free(w->resfac);

  if (w->psi)
    gsl_vector_free(w->psi);

  if (w->dpsi)
    gsl_vector_free(w->dpsi);

  if (w->QSI)
    gsl_matrix_free(w->QSI);

  if (w->D)
    gsl_vector_free(w->D);

  if (w->workn)
    gsl_vector_free(w->workn);

  free(w);
} /* gsl_multifit_robust_free() */

int
gsl_multifit_robust_tune(const double tune, gsl_multifit_robust_workspace *w)
{
  w->tune = tune;
  return GSL_SUCCESS;
}

const char *
gsl_multifit_robust_name(const gsl_multifit_robust_workspace *w)
{
  return w->type->name;
}

gsl_multifit_robust_stats
gsl_multifit_robust_statistics(const gsl_multifit_robust_workspace *w)
{
  return w->stats;
}

/*
gsl_multifit_robust()
  Perform robust iteratively reweighted linear least squares
fit

Inputs: X     - design matrix of basis functions
        y     - right hand side vector
        c     - (output) model coefficients
        cov   - (output) covariance matrix
        w     - workspace
*/

int
gsl_multifit_robust(const gsl_matrix * X,
                    const gsl_vector * y,
                    gsl_vector * c,
                    gsl_matrix * cov,
                    gsl_multifit_robust_workspace *w)
{
  /* check matrix and vector sizes */
  if (X->size1 != y->size)
    {
      GSL_ERROR
        ("number of observations in y does not match rows of matrix X",
         GSL_EBADLEN);
    }
  else if (X->size2 != c->size)
    {
      GSL_ERROR ("number of parameters c does not match columns of matrix X",
                 GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {   
      GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    }   
  else if (c->size != cov->size1)
    {   
      GSL_ERROR
        ("number of parameters does not match size of covariance matrix",
         GSL_EBADLEN);
    }
  else if (X->size1 != w->n || X->size2 != w->p)
    {
      GSL_ERROR
        ("size of workspace does not match size of observation matrix",
         GSL_EBADLEN);
    }
  else
    {
      int s;
      double chisq;
      const double tol = GSL_SQRT_DBL_EPSILON;
      int converged = 0;
      size_t numit = 0;
      const size_t n = y->size;
      double sigy = gsl_stats_sd(y->data, y->stride, n);
      double sig_lower;
      size_t i;

      /*
       * if the initial fit is very good, then finding outliers by comparing
       * them to the residual standard deviation is difficult. Therefore we
       * set a lower bound on the standard deviation estimate that is a small
       * fraction of the standard deviation of the data values
       */
      sig_lower = 1.0e-6 * sigy;
      if (sig_lower == 0.0)
        sig_lower = 1.0;

      /* compute initial estimates using ordinary least squares */
      s = gsl_multifit_linear(X, y, c, cov, &chisq, w->multifit_p);
      if (s)
        return s;

      /* save Q S^{-1} of original matrix */
      gsl_matrix_memcpy(w->QSI, w->multifit_p->QSI);
      gsl_vector_memcpy(w->D, w->multifit_p->D);

      /* compute statistical leverage of each data point */
      s = gsl_linalg_SV_leverage(w->multifit_p->A, w->resfac);
      if (s)
        return s;

      /* correct residuals with factor 1 / sqrt(1 - h) */
      for (i = 0; i < n; ++i)
        {
          double h = gsl_vector_get(w->resfac, i);

          if (h > 0.9999)
            h = 0.9999;

          gsl_vector_set(w->resfac, i, 1.0 / sqrt(1.0 - h));
        }

      /* compute residuals from OLS fit r = y - X c */
      s = gsl_multifit_linear_residuals(X, y, c, w->r);
      if (s)
        return s;

      /* compute estimate of sigma from ordinary least squares */
      w->stats.sigma_ols = gsl_blas_dnrm2(w->r) / sqrt((double) w->stats.dof);

      while (!converged && ++numit <= w->maxiter)
        {
          double sig;

          /* adjust residuals by statistical leverage (see DuMouchel and O'Brien) */
          s = gsl_vector_mul(w->r, w->resfac);
          if (s)
            return s;

          /* compute estimate of standard deviation using MAD */
          sig = robust_madsigma(w->r, w);

          /* scale residuals by standard deviation and tuning parameter */
          gsl_vector_scale(w->r, 1.0 / (GSL_MAX(sig, sig_lower) * w->tune));

          /* compute weights using these residuals */
          s = w->type->wfun(w->r, w->weights);
          if (s)
            return s;

          gsl_vector_memcpy(w->c_prev, c);

          /* solve weighted least squares with new weights */
          s = gsl_multifit_wlinear(X, w->weights, y, c, cov, &chisq, w->multifit_p);
          if (s)
            return s;

          /* compute new residuals r = y - X c */
          s = gsl_multifit_linear_residuals(X, y, c, w->r);
          if (s)
            return s;

          converged = robust_test_convergence(w->c_prev, c, tol);
        }

      /* compute final MAD sigma */
      w->stats.sigma_mad = robust_madsigma(w->r, w);

      /* compute robust estimate of sigma */
      w->stats.sigma_rob = robust_robsigma(w->r, w->stats.sigma_mad, w->tune, w);

      /* compute final estimate of sigma */
      w->stats.sigma = robust_sigma(w->stats.sigma_ols, w->stats.sigma_rob, w);

      /* store number of iterations */
      w->stats.numit = numit;

      {
        double dof = (double) w->stats.dof;
        double rnorm = w->stats.sigma * sqrt(dof); /* see DuMouchel, sec 4.2 */
        double ss_err = rnorm * rnorm;
        double ss_tot = gsl_stats_tss(y->data, y->stride, n);

        /* compute R^2 */
        w->stats.Rsq = 1.0 - ss_err / ss_tot;

        /* compute adjusted R^2 */
        w->stats.adj_Rsq = 1.0 - (1.0 - w->stats.Rsq) * (n - 1.0) / dof;

        /* compute rmse */
        w->stats.rmse = sqrt(ss_err / dof);

        /* store SSE */
        w->stats.sse = ss_err;
      }

      /* calculate covariance matrix = sigma^2 (X^T X)^{-1} */
      s = robust_covariance(w->stats.sigma, cov, w);
      if (s)
        return s;

      /* raise an error if not converged */
      if (numit > w->maxiter)
        {
          GSL_ERROR("maximum iterations exceeded", GSL_EMAXITER);
        }

      return s;
    }
} /* gsl_multifit_robust() */

/* Estimation of values for given x */
int
gsl_multifit_robust_est(const gsl_vector * x, const gsl_vector * c,
                        const gsl_matrix * cov, double *y, double *y_err)
{
  int s = gsl_multifit_linear_est(x, c, cov, y, y_err);

  return s;
}

/***********************************
 * INTERNAL ROUTINES               *
 ***********************************/

/*
robust_test_convergence()
  Test for convergence in robust least squares

Convergence criteria:

|c_i^(k) - c_i^(k-1)| <= tol * max(|c_i^(k)|, |c_i^(k-1)|)

for all i. k refers to iteration number.

Inputs: c_prev - coefficients from previous iteration
        c      - coefficients from current iteration
        tol    - tolerance

Return: 1 if converged, 0 if not
*/

static int
robust_test_convergence(const gsl_vector *c_prev, const gsl_vector *c,
                        const double tol)
{
  size_t p = c->size;
  size_t i;

  for (i = 0; i < p; ++i)
    {
      double ai = gsl_vector_get(c_prev, i);
      double bi = gsl_vector_get(c, i);

      if (fabs(bi - ai) > tol * GSL_MAX(fabs(ai), fabs(bi)))
        return 0; /* not yet converged */
    }

  /* converged */
  return 1;
} /* robust_test_convergence() */

/*
robust_madsigma()
  Estimate the standard deviation of the residuals using
the Median-Absolute-Deviation (MAD) of the residuals,
throwing away the smallest p residuals.

See: Street et al, 1988

Inputs: r - vector of residuals
        w - workspace
*/

static double
robust_madsigma(const gsl_vector *r, gsl_multifit_robust_workspace *w)
{
  gsl_vector_view v;
  double sigma;
  size_t n = r->size;
  const size_t p = w->p;
  size_t i;

  /* copy |r| into workn */
  for (i = 0; i < n; ++i)
    {
      gsl_vector_set(w->workn, i, fabs(gsl_vector_get(r, i)));
    }

  gsl_sort_vector(w->workn);

  /*
   * ignore the smallest p residuals when computing the median
   * (see Street et al 1988)
   */
  v = gsl_vector_subvector(w->workn, p - 1, n - p + 1);
  sigma = gsl_stats_median_from_sorted_data(v.vector.data, v.vector.stride, v.vector.size) / 0.6745;

  return sigma;
} /* robust_madsigma() */

/*
robust_robsigma()
  Compute robust estimate of sigma so that
sigma^2 * inv(X' * X) is a reasonable estimate of
the covariance for robust regression. Based heavily
on the equations of Street et al, 1988.

Inputs: r    - vector of residuals y - X c
        s    - sigma estimate using MAD
        tune - tuning constant
        w    - workspace
*/

static double
robust_robsigma(const gsl_vector *r, const double s,
                const double tune, gsl_multifit_robust_workspace *w)
{
  double sigma;
  size_t i;
  const size_t n = w->n;
  const size_t p = w->p;
  const double st = s * tune;
  double a, b, lambda;

  /* compute u = r / sqrt(1 - h) / st */
  gsl_vector_memcpy(w->workn, r);
  gsl_vector_mul(w->workn, w->resfac);
  gsl_vector_scale(w->workn, 1.0 / st);

  /* compute w(u) and psi'(u) */
  w->type->wfun(w->workn, w->psi);
  w->type->psi_deriv(w->workn, w->dpsi);

  /* compute psi(u) = u*w(u) */
  gsl_vector_mul(w->psi, w->workn);

  /* Street et al, Eq (3) */
  a = gsl_stats_mean(w->dpsi->data, w->dpsi->stride, n);

  /* Street et al, Eq (5) */
  b = 0.0;
  for (i = 0; i < n; ++i)
    {
      double psi_i = gsl_vector_get(w->psi, i);
      double resfac = gsl_vector_get(w->resfac, i);
      double fac = 1.0 / (resfac*resfac); /* 1 - h */

      b += fac * psi_i * psi_i;
    }
  b /= (double) (n - p);

  /* Street et al, Eq (5) */
  lambda = 1.0 + ((double)p)/((double)n) * (1.0 - a) / a;

  sigma = lambda * sqrt(b) * st / a;

  return sigma;
} /* robust_robsigma() */

/*
robust_sigma()
  Compute final estimate of residual standard deviation, using
the OLS and robust sigma estimates.

This equation is taken from DuMouchel and O'Brien, sec 4.1:
\hat{\sigma_R}

Inputs: s_ols - OLS sigma
        s_rob - robust sigma
        w     - workspace

Return: final estimate of sigma
*/

static double
robust_sigma(const double s_ols, const double s_rob,
             gsl_multifit_robust_workspace *w)
{
  double sigma;
  const size_t p = w->p;
  const size_t n = w->n;

  /* see DuMouchel and O'Brien, sec 4.1 */
  sigma = GSL_MAX(s_rob,
                  sqrt((s_ols*s_ols*p*p + s_rob*s_rob*n) /
                       (p*p + n)));

  return sigma;
} /* robust_sigma() */

/*
robust_covariance()
  Calculate final covariance matrix, defined as:

  sigma * (X^T X)^{-1}

Inputs: sigma - residual standard deviation
        cov   - (output) covariance matrix
        w     - workspace
*/

static int
robust_covariance(const double sigma, gsl_matrix *cov,
                  gsl_multifit_robust_workspace *w)
{
  int s = 0;
  const size_t p = w->p;
  const double s2 = sigma * sigma;
  size_t i, j;
  gsl_matrix *QSI = w->QSI;
  gsl_vector *D = w->D;

  /* Form variance-covariance matrix cov = s2 * (Q S^-1) (Q S^-1)^T */

  for (i = 0; i < p; i++)
    {
      gsl_vector_view row_i = gsl_matrix_row (QSI, i);
      double d_i = gsl_vector_get (D, i);

      for (j = i; j < p; j++)
        {
          gsl_vector_view row_j = gsl_matrix_row (QSI, j);
          double d_j = gsl_vector_get (D, j);
          double s;

          gsl_blas_ddot (&row_i.vector, &row_j.vector, &s);

          gsl_matrix_set (cov, i, j, s * s2 / (d_i * d_j));
          gsl_matrix_set (cov, j, i, s * s2 / (d_i * d_j));
        }
    }

  return s;
} /* robust_covariance() */
