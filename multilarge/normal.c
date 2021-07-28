/* normal.c
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>

typedef struct
{
  size_t p;              /* number of columns of LS matrix */
  gsl_matrix *ATA;       /* A^T A, p-by-p */
  gsl_vector *ATb;       /* A^T b, p-by-1 */
  double normb;          /* || b || */
  gsl_matrix *work_ATA;  /* workspace for chol(ATA), p-by-p */
  gsl_vector *workp;     /* workspace size p */
  gsl_vector *work3p;    /* workspace size 3*p */
  gsl_vector *D;         /* scale factors for ATA, size p */
  gsl_vector *c;         /* solution vector for L-curve */
  int eigen;             /* 1 if eigenvalues computed */
  double eval_min;       /* minimum eigenvalue */
  double eval_max;       /* maximum eigenvalue */
  gsl_eigen_symm_workspace *eigen_p;
} normal_state_t;

static void *normal_alloc(const size_t p);
static void normal_free(void *vstate);
static int normal_reset(void *vstate);
static int normal_accumulate(gsl_matrix * A, gsl_vector * b,
                             void * vstate);
static int normal_solve(const double lambda, gsl_vector * x,
                        double * rnorm, double * snorm,
                        void * vstate);
static int normal_rcond(double * rcond, void * vstate);
static int normal_lcurve(gsl_vector * reg_param, gsl_vector * rho,
                         gsl_vector * eta, void * vstate);
static const gsl_matrix * normal_ATA(const void * vstate);
static const gsl_vector * normal_ATb(const void * vstate);
static int normal_solve_system(const double lambda, gsl_vector * x,
                               normal_state_t *state);
static int normal_solve_cholesky(gsl_matrix * ATA, const gsl_vector * ATb,
                                 gsl_vector * x, normal_state_t *state);
static int normal_calc_norms(const gsl_vector *x, double *rnorm,
                             double *snorm, normal_state_t *state);
static int normal_eigen(normal_state_t *state);

/*
normal_alloc()
  Allocate workspace for solving large linear least squares
problems using the normal equations approach

Inputs: p    - number of columns of LS matrix

Return: pointer to workspace
*/

static void *
normal_alloc(const size_t p)
{
  normal_state_t *state;

  if (p == 0)
    {
      GSL_ERROR_NULL("p must be a positive integer",
                     GSL_EINVAL);
    }

  state = calloc(1, sizeof(normal_state_t));
  if (!state)
    {
      GSL_ERROR_NULL("failed to allocate normal state", GSL_ENOMEM);
    }

  state->p = p;

  state->ATA = gsl_matrix_alloc(p, p);
  if (state->ATA == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate ATA matrix", GSL_ENOMEM);
    }

  state->work_ATA = gsl_matrix_alloc(p, p);
  if (state->work_ATA == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate temporary ATA matrix", GSL_ENOMEM);
    }

  state->ATb = gsl_vector_alloc(p);
  if (state->ATb == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate ATb vector", GSL_ENOMEM);
    }

  state->D = gsl_vector_alloc(p);
  if (state->D == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate D vector", GSL_ENOMEM);
    }

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate temporary ATb vector", GSL_ENOMEM);
    }

  state->work3p = gsl_vector_alloc(3 * p);
  if (state->work3p == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate work3p", GSL_ENOMEM);
    }

  state->c = gsl_vector_alloc(p);
  if (state->c == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate c vector", GSL_ENOMEM);
    }

  state->eigen_p = gsl_eigen_symm_alloc(p);
  if (state->eigen_p == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate eigen workspace", GSL_ENOMEM);
    }

  normal_reset(state);

  return state;
}

static void
normal_free(void *vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  if (state->ATA)
    gsl_matrix_free(state->ATA);

  if (state->work_ATA)
    gsl_matrix_free(state->work_ATA);

  if (state->ATb)
    gsl_vector_free(state->ATb);

  if (state->D)
    gsl_vector_free(state->D);

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->work3p)
    gsl_vector_free(state->work3p);

  if (state->c)
    gsl_vector_free(state->c);

  if (state->eigen_p)
    gsl_eigen_symm_free(state->eigen_p);

  free(state);
}

static int
normal_reset(void *vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  gsl_matrix_set_zero(state->ATA);
  gsl_vector_set_zero(state->ATb);
  state->normb = 0.0;
  state->eigen = 0;
  state->eval_min = 0.0;
  state->eval_max = 0.0;

  return GSL_SUCCESS;
}

/*
normal_accumulate()
  Add a new block of rows to the normal equations system

Inputs: A      - new block of rows, n-by-p
        b      - new rhs vector n-by-1
        vstate - workspace

Return: success/error
*/

static int
normal_accumulate(gsl_matrix * A, gsl_vector * b, void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;
  const size_t n = A->size1;

  if (A->size2 != state->p)
    {
      GSL_ERROR("columns of A do not match workspace", GSL_EBADLEN);
    }
  else if (n != b->size)
    {
      GSL_ERROR("A and b have different numbers of rows", GSL_EBADLEN);
    }
  else
    {
      int s;

      /* ATA += A^T A, using only the lower half of the matrix */
      s = gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, A, 1.0, state->ATA);
      if (s)
        return s;

      /* ATb += A^T b */
      s = gsl_blas_dgemv(CblasTrans, 1.0, A, b, 1.0, state->ATb);
      if (s)
        return s;

      /* update || b || */
      state->normb = gsl_hypot(state->normb, gsl_blas_dnrm2(b));

      return GSL_SUCCESS;
    }
}

/*
normal_solve()
  Solve normal equations system:

(A^T A + \lambda^2 I) x = A^T b

using Cholesky decomposition

Inputs: lambda - regularization parameter
        x      - (output) solution vector p-by-1
        rnorm  - (output) residual norm ||b - A x||
        snorm  - (output) solution norm ||x||
        vstate - workspace

Return: success/error
*/

static int
normal_solve(const double lambda, gsl_vector * x,
             double * rnorm, double * snorm,
             void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  if (x->size != state->p)
    {
      GSL_ERROR("solution vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* solve system (A^T A) x = A^T b */
      status = normal_solve_system(lambda, x, state);
      if (status)
        {
          GSL_ERROR("failed to solve normal equations", status);
        }

      /* compute residual norm ||y - X c|| and solution norm ||x|| */
      normal_calc_norms(x, rnorm, snorm, state);

      return GSL_SUCCESS;
    }
}

static int
normal_rcond(double * rcond, void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;
  int status = GSL_SUCCESS;
  double rcond_ATA;

  status = gsl_linalg_cholesky_rcond(state->work_ATA, &rcond_ATA, state->work3p);
  if (status == GSL_SUCCESS)
    *rcond = sqrt(rcond_ATA);

  return status;
}

/*
normal_lcurve()
  Compute L-curve of least squares system

Inputs: reg_param - (output) vector of regularization parameters
        rho       - (output) vector of residual norms
        eta       - (output) vector of solution norms
        vstate    - workspace

Return: success/error
*/

static int
normal_lcurve(gsl_vector * reg_param, gsl_vector * rho,
              gsl_vector * eta, void * vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;
  int status;
  double smin, smax; /* minimum/maximum singular values */
  size_t i;

  if (state->eigen == 0)
    {
      status = normal_eigen(state);
      if (status)
        return status;
    }

  if (state->eval_max < 0.0)
    {
      GSL_ERROR("matrix is not positive definite", GSL_EDOM);
    }

  /* compute singular values which are sqrts of eigenvalues */
  smax = sqrt(state->eval_max);
  if (state->eval_min > 0.0)
    smin = sqrt(state->eval_min);
  else
    smin = 0.0;

  /* compute vector of regularization parameters */
  gsl_multifit_linear_lreg(smin, smax, reg_param);

  /* solve normal equations for each regularization parameter */
  for (i = 0; i < reg_param->size; ++i)
    {
      double lambda = gsl_vector_get(reg_param, i);
      double rnorm, snorm;

      status = normal_solve_system(lambda, state->c, state);
      if (status)
        return status;

      /* compute ||y - X c|| and ||c|| */
      normal_calc_norms(state->c, &rnorm, &snorm, state);

      gsl_vector_set(rho, i, rnorm);
      gsl_vector_set(eta, i, snorm);
    }

  return GSL_SUCCESS;
}

static const gsl_matrix *
normal_ATA(const void * vstate)
{
  const normal_state_t *state = (const normal_state_t *) vstate;
  return state->ATA;
}

static const gsl_vector *
normal_ATb(const void * vstate)
{
  const normal_state_t *state = (const normal_state_t *) vstate;
  return state->ATb;
}

/*
normal_solve_system()
  Compute solution to normal equations:

(A^T A + lambda^2*I) x = A^T b

using LDL decomposition.

Inputs: x     - (output) solution vector
        state - workspace

Return: success/error
*/

static int
normal_solve_system(const double lambda, gsl_vector * x, normal_state_t *state)
{
  int status;
  const double lambda_sq = lambda * lambda;
  gsl_vector_view d = gsl_matrix_diagonal(state->work_ATA);

  /* copy ATA matrix to temporary workspace and regularize */
  gsl_matrix_tricpy(CblasLower, CblasNonUnit, state->work_ATA, state->ATA);
  gsl_vector_add_constant(&d.vector, lambda_sq);

  /* solve with Cholesky decomposition */
  status = normal_solve_cholesky(state->work_ATA, state->ATb, x, state);
  if (status)
    return status;

  return status;
}

static int
normal_solve_cholesky(gsl_matrix * ATA, const gsl_vector * ATb,
                      gsl_vector * x, normal_state_t *state)
{
  int status;

  status = gsl_linalg_cholesky_decomp2(ATA, state->D);
  if (status)
    return status;

  status = gsl_linalg_cholesky_solve2(ATA, state->D, ATb, x);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
normal_calc_norms()
  Compute residual norm ||y - X c|| and solution
norm ||c||

Inputs: x     - solution vector
        rnorm - (output) residual norm ||y - X c||
        snorm - (output) solution norm ||c||
        state - workspace
*/

static int
normal_calc_norms(const gsl_vector *x, double *rnorm,
                  double *snorm, normal_state_t *state)
{
  double r2;

  /* compute solution norm ||x|| */
  *snorm = gsl_blas_dnrm2(x);

  /* compute residual norm ||b - Ax|| */

  /* compute: A^T A x - 2 A^T b */
  gsl_vector_memcpy(state->workp, state->ATb);
  gsl_blas_dsymv(CblasLower, 1.0, state->ATA, x, -2.0, state->workp);

  /* compute: x^T A^T A x - 2 x^T A^T b */
  gsl_blas_ddot(x, state->workp, &r2);

  /* add b^T b */
  r2 += state->normb * state->normb;

  *rnorm = sqrt(r2);

  return GSL_SUCCESS;
}

/*
normal_eigen()
  Compute eigenvalues of A^T A matrix, which
are stored in state->workp on output. Also,
state->eval_min and state->eval_max are set
to the minimum/maximum eigenvalues
*/

static int
normal_eigen(normal_state_t *state)
{
  int status;

  /* copy lower triangle of ATA to temporary workspace */
  gsl_matrix_tricpy(CblasLower, CblasNonUnit, state->work_ATA, state->ATA);

  /* compute eigenvalues of ATA */
  status = gsl_eigen_symm(state->work_ATA, state->workp, state->eigen_p);
  if (status)
    return status;

  gsl_vector_minmax(state->workp, &state->eval_min, &state->eval_max);

  state->eigen = 1;

  return GSL_SUCCESS;
}

static const gsl_multilarge_linear_type normal_type =
{
  "normal",
  normal_alloc,
  normal_reset,
  normal_accumulate,
  normal_solve,
  normal_rcond,
  normal_lcurve,
  normal_ATA,
  normal_ATb,
  normal_free
};

const gsl_multilarge_linear_type * gsl_multilarge_linear_normal =
  &normal_type;
