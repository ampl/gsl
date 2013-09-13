/* multifit/gsl_multifit.h
 * 
 * Copyright (C) 2000, 2007, 2010 Brian Gough
 * Copyright (C) 2013, Patrick Alken
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

#ifndef __GSL_MULTIFIT_H__
#define __GSL_MULTIFIT_H__

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_types.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct 
{
  size_t n; /* number of observations */
  size_t p; /* number of parameters */
  gsl_matrix * A;
  gsl_matrix * Q;
  gsl_matrix * QSI;
  gsl_vector * S;
  gsl_vector * t;
  gsl_vector * xt;
  gsl_vector * D;
} 
gsl_multifit_linear_workspace;

gsl_multifit_linear_workspace *
gsl_multifit_linear_alloc (size_t n, size_t p);

void
gsl_multifit_linear_free (gsl_multifit_linear_workspace * work);

int
gsl_multifit_linear (const gsl_matrix * X,
                     const gsl_vector * y,
                     gsl_vector * c,
                     gsl_matrix * cov,
                     double * chisq,
                     gsl_multifit_linear_workspace * work);

int
gsl_multifit_linear_svd (const gsl_matrix * X,
                         const gsl_vector * y,
                         double tol,
                         size_t * rank,
                         gsl_vector * c,
                         gsl_matrix * cov,
                         double *chisq, 
                         gsl_multifit_linear_workspace * work);

int
gsl_multifit_linear_usvd (const gsl_matrix * X,
                          const gsl_vector * y,
                          double tol,
                          size_t * rank,
                          gsl_vector * c,
                          gsl_matrix * cov,
                          double *chisq, 
                          gsl_multifit_linear_workspace * work);

int
gsl_multifit_wlinear (const gsl_matrix * X,
                      const gsl_vector * w,
                      const gsl_vector * y,
                      gsl_vector * c,
                      gsl_matrix * cov,
                      double * chisq,
                      gsl_multifit_linear_workspace * work);

int
gsl_multifit_wlinear_svd (const gsl_matrix * X,
                          const gsl_vector * w,
                          const gsl_vector * y,
                          double tol,
                          size_t * rank,
                          gsl_vector * c,
                          gsl_matrix * cov,
                          double *chisq, 
                          gsl_multifit_linear_workspace * work);

int
gsl_multifit_wlinear_usvd (const gsl_matrix * X,
                           const gsl_vector * w,
                           const gsl_vector * y,
                           double tol,
                           size_t * rank,
                           gsl_vector * c,
                           gsl_matrix * cov,
                           double *chisq, 
                           gsl_multifit_linear_workspace * work);

int
gsl_multifit_linear_est (const gsl_vector * x,
                         const gsl_vector * c,
                         const gsl_matrix * cov, double *y, double *y_err);

int
gsl_multifit_linear_residuals (const gsl_matrix *X, const gsl_vector *y,
                               const gsl_vector *c, gsl_vector *r);

typedef struct
{
  const char * name;     /* method name */
  int (*wfun)(const gsl_vector *r, gsl_vector *w);
  int (*psi_deriv)(const gsl_vector *r, gsl_vector *dpsi);
  double tuning_default; /* default tuning constant */
} gsl_multifit_robust_type;

typedef struct
{
  double sigma_ols;    /* OLS estimate of sigma */
  double sigma_mad;    /* MAD estimate of sigma */
  double sigma_rob;    /* robust estimate of sigma */
  double sigma;        /* final estimate of sigma */
  double Rsq;          /* R^2 coefficient of determination */
  double adj_Rsq;      /* degree of freedom adjusted R^2 */
  double rmse;         /* root mean squared error */
  double sse;          /* residual sum of squares */
  size_t dof;          /* degrees of freedom */
  size_t numit;        /* number of iterations */
  gsl_vector *weights; /* final weights */
  gsl_vector *r;       /* final residuals y - X c */
} gsl_multifit_robust_stats;

typedef struct
{
  size_t n;            /* number of observations */
  size_t p;            /* number of parameters */
  size_t numit;        /* number of iterations */
  size_t maxiter;      /* maximum iterations */
  const gsl_multifit_robust_type *type;
  double tune;         /* tuning parameter */

  gsl_vector *r;       /* residuals at current iteration */
  gsl_vector *weights; /* weights at current iteration */
  gsl_vector *c_prev;  /* coefficients from previous iteration */
  gsl_vector *resfac;  /* multiplicative factors for residuals */

  gsl_vector *psi;     /* psi(r) */
  gsl_vector *dpsi;    /* psi'(r) */

  gsl_matrix *QSI;     /* Q S^{-1} of original matrix X */
  gsl_vector *D;       /* balancing parameters of original matrix X */

  gsl_vector *workn;   /* workspace of length n */

  gsl_multifit_robust_stats stats; /* various statistics */

  gsl_multifit_linear_workspace *multifit_p;
} gsl_multifit_robust_workspace;

/* available types */
GSL_VAR const gsl_multifit_robust_type * gsl_multifit_robust_default;
GSL_VAR const gsl_multifit_robust_type * gsl_multifit_robust_bisquare;
GSL_VAR const gsl_multifit_robust_type * gsl_multifit_robust_cauchy;
GSL_VAR const gsl_multifit_robust_type * gsl_multifit_robust_fair;
GSL_VAR const gsl_multifit_robust_type * gsl_multifit_robust_huber;
GSL_VAR const gsl_multifit_robust_type * gsl_multifit_robust_ols;
GSL_VAR const gsl_multifit_robust_type * gsl_multifit_robust_welsch;

gsl_multifit_robust_workspace *gsl_multifit_robust_alloc(const gsl_multifit_robust_type *T,
                                                         const size_t n, const size_t p);
void gsl_multifit_robust_free(gsl_multifit_robust_workspace *w);
int gsl_multifit_robust_tune(const double tune, gsl_multifit_robust_workspace *w);
const char *gsl_multifit_robust_name(const gsl_multifit_robust_workspace *w);
gsl_multifit_robust_stats gsl_multifit_robust_statistics(const gsl_multifit_robust_workspace *w);
int gsl_multifit_robust(const gsl_matrix * X, const gsl_vector * y,
                        gsl_vector * c, gsl_matrix *cov,
                        gsl_multifit_robust_workspace *w);
int gsl_multifit_robust_est(const gsl_vector * x, const gsl_vector * c,
                            const gsl_matrix * cov, double *y, double *y_err);

__END_DECLS

#endif /* __GSL_MULTIFIT_H__ */
