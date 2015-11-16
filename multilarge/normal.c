/* normal.c
 * 
 * Copyright (C) 2015 Patrick Alken
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
#include <gsl/gsl_multilarge.h>

typedef struct
{
  size_t p;             /* number of columns of LS matrix */
  gsl_matrix *ATA;      /* A^T A, p-by-p */
  gsl_vector *ATb;      /* A^T b, p-by-1 */
  double normb;         /* || b || */
  gsl_matrix *work_ATA; /* workspace for chol(ATA), p-by-p */
  gsl_vector *workp;    /* workspace size p */
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
static int normal_solve_system(const double lambda, gsl_vector * x,
                               normal_state_t *state);
static int normal_solve_cholesky(gsl_matrix * ATA, const gsl_vector * ATb,
                                 gsl_vector * x, normal_state_t *state);
static int normal_solve_QR(gsl_matrix * ATA, const gsl_vector * ATb,
                           gsl_vector * x, normal_state_t *state);
static int normal_copy_lower(gsl_matrix * dest, const gsl_matrix * src);
static int normal_copy_lowup(gsl_matrix * A);

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

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      normal_free(state);
      GSL_ERROR_NULL("failed to allocate temporary ATb vector", GSL_ENOMEM);
    }

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

  if (state->workp)
    gsl_vector_free(state->workp);

  free(state);
}

static int
normal_reset(void *vstate)
{
  normal_state_t *state = (normal_state_t *) vstate;

  gsl_matrix_set_zero(state->ATA);
  gsl_vector_set_zero(state->ATb);
  state->normb = 0.0;

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
      double r2;

      /* solve system (A^T A) x = A^T b */
      status = normal_solve_system(lambda, x, state);
      if (status)
        {
          GSL_ERROR("failed to solve normal equations", status);
        }

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
}

static int
normal_rcond(double * rcond, void * vstate)
{
  *rcond = 0.0;
  return GSL_SUCCESS;
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
  return GSL_SUCCESS;
}

/*
normal_solve_system()
  Compute solution to normal equations:

(A^T A + lambda^2*I) x = A^T b

First we try Cholesky decomposition. If that
fails, try QR

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
  gsl_error_handler_t *err_handler;

  /* copy ATA matrix to temporary workspace and regularize */
  normal_copy_lower(state->work_ATA, state->ATA);
  gsl_vector_add_constant(&d.vector, lambda_sq);

  /* turn off error handler in case Cholesky fails */
  err_handler = gsl_set_error_handler_off();

  /* try Cholesky decomposition first */
  status = normal_solve_cholesky(state->work_ATA, state->ATb, x, state);
  if (status)
    {
      /* restore ATA matrix and try QR decomposition */
      normal_copy_lower(state->work_ATA, state->ATA);
      gsl_vector_add_constant(&d.vector, lambda_sq);

      status = normal_solve_QR(state->work_ATA, state->ATb, x, state);
    }

  /* restore error handler */
  gsl_set_error_handler(err_handler);

  return status;
}

static int
normal_solve_cholesky(gsl_matrix * ATA, const gsl_vector * ATb,
                      gsl_vector * x, normal_state_t *state)
{
  int status;

  /* compute Cholesky decomposition of A^T A */
  status = gsl_linalg_cholesky_decomp(ATA);
  if (status)
    return status;

  status = gsl_linalg_cholesky_solve(ATA, ATb, x);
  if (status)
    return status;

  return GSL_SUCCESS;
}

static int
normal_solve_QR(gsl_matrix * ATA, const gsl_vector * ATb,
                gsl_vector * x, normal_state_t *state)
{
  int status;

  /* copy lower triangle of ATA to upper */
  normal_copy_lowup(ATA);

  status = gsl_linalg_QR_decomp(ATA, state->workp);
  if (status)
    return status;

  status = gsl_linalg_QR_solve(ATA, state->workp, ATb, x);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/* copy lower triangle of src to dest, including diagonal */
static int
normal_copy_lower(gsl_matrix * dest, const gsl_matrix * src)
{
  const size_t src_size1 = src->size1;
  const size_t src_tda = src->tda;
  const size_t dest_tda = dest->tda;
  size_t i, j;

  for (i = 0; i < src_size1 ; i++)
    {
      for (j = 0; j <= i; j++)
        {
          dest->data[dest_tda * i + j] 
            = src->data[src_tda * i + j];
        }
    }

  return GSL_SUCCESS;
}

/* copy lower triangular part of matrix to upper */
static int
normal_copy_lowup(gsl_matrix * A)
{
  const size_t size1 = A->size1;
  const size_t size2 = A->size2;
  const size_t tda = A->tda;
  size_t i, j;

  for (i = 0; i < size1; ++i)
    {
      for (j = i + 1; j < size2; ++j)
        {
          A->data[tda * i + j] 
            = A->data[tda * j + i];
        }
    }

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
  normal_free
};

const gsl_multilarge_linear_type * gsl_multilarge_linear_normal =
  &normal_type;
