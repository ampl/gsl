/* tsqr.c
 * 
 * Copyright (C) 2015, 2019 Patrick Alken
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
 * This module implements the sequential TSQR algorithm
 * described in
 *
 * [1] Demmel, J., Grigori, L., Hoemmen, M. F., and Langou, J.
 *     "Communication-optimal parallel and sequential QR and LU factorizations",
 *     UCB Technical Report No. UCB/EECS-2008-89, 2008.
 *
 * The algorithm operates on a tall least squares system:
 *
 * [ A_1 ] x = [ b_1 ]
 * [ A_2 ]     [ b_2 ]
 * [ ... ]     [ ... ]
 * [ A_k ]     [ b_k ]
 *
 * as follows:
 *
 * 1. Initialize
 *    a. [Q_1,R_1] = qr(A_1)
 *    b. y_1 = Q_1^T b_1
 * 2. Loop i = 2:k
 *    a. [Q_i,R_i] = qr( [ R_{i-1} ; A_i ] )
 *    b. y_i = Q_i^T [ y_{i-1} ; b_i ]
 * 3. Output:
 *    a. R = R_k
 *    b. Q^T b = y_k
 *
 * Step 2(a) is optimized to take advantage
 * of the sparse structure of the matrix.
 *
 * We can write the full n-by-1 vector,
 *
 * Q^T b = [ z_1 ] p
 *         [ z_2 ] n - p
 *
 * The first p elements z_1 are stored on each accumulation,
 * while only the norm ||z_2|| is stored on each accumulation.
 * The norm ||z_2|| is required to calculate the final residual norm.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>

typedef struct
{
  size_t p;             /* number of columns of LS matrix */
  int nblocks;          /* number of blocks processed */
  int svd;              /* SVD of R factor has been computed */

  gsl_matrix *T;        /* block reflector matrix, p-by-p */
  gsl_matrix *R;        /* R factor, p-by-p */
  gsl_vector *QTb;      /* [ (Q^T b)(1:p) ; ||z_2|| ], size (p+1)-by-1 */
  gsl_vector *work;     /* workspace, size p */
  gsl_vector *work3;    /* workspace, size 3*p */

  gsl_multifit_linear_workspace * multifit_workspace_p;
} tsqr_state_t;

static void *tsqr_alloc(const size_t p);
static void tsqr_free(void *vstate);
static int tsqr_reset(void *vstate);
static int tsqr_accumulate(gsl_matrix * A, gsl_vector * b,
                           void * vstate);
static int tsqr_solve(const double lambda, gsl_vector * x,
                      double * rnorm, double * snorm,
                      void * vstate);
static int tsqr_rcond(double * rcond, void * vstate);
static int tsqr_lcurve(gsl_vector * reg_param, gsl_vector * rho,
                       gsl_vector * eta, void * vstate);
static const gsl_matrix * tsqr_R(const void * vstate);
static const gsl_vector * tsqr_QTb(const void * vstate);
static int tsqr_svd(tsqr_state_t * state);

/*
tsqr_alloc()
  Allocate workspace for solving large linear least squares
problems using the TSQR approach

Inputs: p - number of columns of LS matrix

Return: pointer to workspace
*/

static void *
tsqr_alloc(const size_t p)
{
  tsqr_state_t *state;

  if (p == 0)
    {
      GSL_ERROR_NULL("p must be a positive integer",
                     GSL_EINVAL);
    }

  state = calloc(1, sizeof(tsqr_state_t));
  if (!state)
    {
      GSL_ERROR_NULL("failed to allocate tsqr state", GSL_ENOMEM);
    }

  state->p = p;
  state->nblocks = 0;

  state->R = gsl_matrix_alloc(p, p);
  if (state->R == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate R matrix", GSL_ENOMEM);
    }

  state->QTb = gsl_vector_calloc(p + 1);
  if (state->QTb == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate QTb vector", GSL_ENOMEM);
    }

  state->T = gsl_matrix_alloc(p, p);
  if (state->T == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate T matrix", GSL_ENOMEM);
    }

  state->work = gsl_vector_alloc(p);
  if (state->work == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate workspace vector", GSL_ENOMEM);
    }

  state->work3 = gsl_vector_alloc(3 * p);
  if (state->work3 == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate work3 vector", GSL_ENOMEM);
    }

  state->multifit_workspace_p = gsl_multifit_linear_alloc(p, p);
  if (state->multifit_workspace_p == NULL)
    {
      tsqr_free(state);
      GSL_ERROR_NULL("failed to allocate multifit workspace", GSL_ENOMEM);
    }

  return state;
}

static void
tsqr_free(void *vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;

  if (state->R)
    gsl_matrix_free(state->R);

  if (state->QTb)
    gsl_vector_free(state->QTb);

  if (state->T)
    gsl_matrix_free(state->T);

  if (state->work)
    gsl_vector_free(state->work);

  if (state->work3)
    gsl_vector_free(state->work3);

  if (state->multifit_workspace_p)
    gsl_multifit_linear_free(state->multifit_workspace_p);

  free(state);
}

static int
tsqr_reset(void *vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;

  gsl_matrix_set_zero(state->R);
  gsl_vector_set_zero(state->QTb);
  state->nblocks = 0;
  state->svd = 0;

  return GSL_SUCCESS;
}

/*
tsqr_accumulate()
  Add a new block of rows to the QR system

Inputs: A      - new block of rows, n-by-p
        b      - new rhs vector n-by-1
        vstate - workspace

Return: success/error

Notes:
1) On output, the upper triangular portion of state->R(1:p,1:p)
contains current R matrix

2) state->QTb(1:p) contains the first p elements of Q^T b

3) state->QTb(p+1) contains ||z_2||

4) A and b are destroyed
*/

static int
tsqr_accumulate(gsl_matrix * A, gsl_vector * b, void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;
  const size_t n = A->size1;
  const size_t p = A->size2;

  if (p != state->p)
    {
      GSL_ERROR("columns of A do not match workspace", GSL_EBADLEN);
    }
  else if (n != b->size)
    {
      GSL_ERROR("A and b have different numbers of rows", GSL_EBADLEN);
    }
  else if (state->nblocks == 0 && n < p)
    {
      GSL_ERROR ("n must be >= p", GSL_EBADLEN);
    }
  else if (state->nblocks == 0)
    {
      int status;
      gsl_matrix_view R = gsl_matrix_submatrix(A, 0, 0, p, p);
      gsl_vector_view QTb0 = gsl_vector_subvector(state->QTb, 0, p);
      double * norm_z2 = gsl_vector_ptr(state->QTb, p);
      gsl_vector_view b1 = gsl_vector_subvector(b, 0, p);

      /* this is the first matrix block A_1, compute its (dense) QR decomposition */

      /* compute QR decomposition of A */
      status = gsl_linalg_QR_decomp_r(A, state->T);
      if (status)
        return status;

      /* store upper triangular R factor in state->R */
      gsl_matrix_tricpy(CblasUpper, CblasNonUnit, state->R, &R.matrix);

      /* compute Q^T b and keep the first p elements */
      gsl_linalg_QR_QTvec_r(A, state->T, b, state->work);
      gsl_vector_memcpy(&QTb0.vector, &b1.vector);

      if (n > p)
        {
          gsl_vector_view b2 = gsl_vector_subvector(b, p, n - p);
          *norm_z2 = gsl_blas_dnrm2(&b2.vector);
        }
      else
        *norm_z2 = 0.0;

      state->nblocks = 1;

      return GSL_SUCCESS;
    }
  else
    {
      int status;
      gsl_vector_view QTb0 = gsl_vector_subvector(state->QTb, 0, p);
      double * norm_z2 = gsl_vector_ptr(state->QTb, p);

      /* compute QR decomposition of [ R_{i-1} ; A_i ], accounting for
       * sparse structure */
      status = gsl_linalg_QR_UR_decomp(state->R, A, state->T);
      if (status)
        return status;

      /*
       * Compute:
       *
       * Q^T [ QTb_{i-1} ] = [ QTb_{i-1} - w ]
       *     [    b_i    ]   [  b_i - V~ w   ]
       *
       * where:
       * 
       * w = T^T (QTb_{i-1} + V~^T b_i)
       *
       *       p
       * V = [ I  ] p
       *     [ V~ ] n
       */
      gsl_vector_memcpy(state->work, &QTb0.vector);
      gsl_blas_dgemv(CblasTrans, 1.0, A, b, 1.0, state->work);                     /* w := w + V~^T b */
      gsl_blas_dtrmv(CblasUpper, CblasTrans, CblasNonUnit, state->T, state->work); /* w := T^T w */
      gsl_vector_sub(&QTb0.vector, state->work);                                   /* QTb := QTb - w */

      /* update residual norm */
      gsl_blas_dgemv(CblasNoTrans, -1.0, A, state->work, 1.0, b);                  /* b := b - V~ w */
      *norm_z2 = gsl_hypot(*norm_z2, gsl_blas_dnrm2(b));

      return GSL_SUCCESS;
    }
}

/*
tsqr_solve()
  Solve the least squares system:

chi^2 = || QTb - R x ||^2 + lambda^2 || x ||^2

using the SVD of R

Inputs: lambda - regularization parameter
        x      - (output) solution vector p-by-1
        rnorm  - (output) residual norm ||b - A x||
        snorm  - (output) solution norm ||x||
        vstate - workspace

Return: success/error
*/

static int
tsqr_solve(const double lambda, gsl_vector * x,
           double * rnorm, double * snorm,
           void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;

  if (x->size != state->p)
    {
      GSL_ERROR ("solution vector does not match workspace", GSL_EBADLEN);
    }
  else if (lambda < 0.0)
    {
      GSL_ERROR ("regularization parameter should be non-negative", GSL_EINVAL);
    }
  else
    {
      gsl_vector_view QTb0 = gsl_vector_subvector(state->QTb, 0, state->p);
      double norm_z2 = gsl_vector_get(state->QTb, state->p); /* || z_2 || */

      if (lambda == 0.0)
        {
          /* solve: R x = Q^T b */
          gsl_vector_memcpy(x, &QTb0.vector);
          gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, state->R, x);
          *rnorm = norm_z2;
          *snorm = gsl_blas_dnrm2(x);
        }
      else
        {
          int status;

          /* compute SVD of R if not already computed */
          if (state->svd == 0)
            {
              status = tsqr_svd(state);
              if (status)
                return status;
            }

          status = gsl_multifit_linear_solve(lambda, state->R, &QTb0.vector, x, rnorm, snorm,
                                             state->multifit_workspace_p);
          if (status)
            return status;

          *rnorm = gsl_hypot(*rnorm, norm_z2);
        }

      return GSL_SUCCESS;
    }
}

/*
tsqr_lcurve()
  Compute L-curve of least squares system

Inputs: reg_param - (output) vector of regularization parameters
        rho       - (output) vector of residual norms
        eta       - (output) vector of solution norms
        vstate    - workspace

Return: success/error
*/

static int
tsqr_lcurve(gsl_vector * reg_param, gsl_vector * rho,
            gsl_vector * eta, void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;
  gsl_vector_view QTb0 = gsl_vector_subvector(state->QTb, 0, state->p);
  double norm_z2 = gsl_vector_get(state->QTb, state->p);
  int status;
  size_t i;

  /* compute SVD of R if not already computed */
  if (state->svd == 0)
    {
      status = tsqr_svd(state);
      if (status)
        return status;
    }

  status = gsl_multifit_linear_lcurve(&QTb0.vector, reg_param, rho, eta,
                                      state->multifit_workspace_p);

  /* now add contribution to rnorm from Q2 factor ||z_2|| */
  for (i = 0; i < rho->size; ++i)
    {
      double *rhoi = gsl_vector_ptr(rho, i);
      *rhoi = gsl_hypot(*rhoi, norm_z2);
    }

  return status;
}

static const gsl_matrix *
tsqr_R(const void * vstate)
{
  const tsqr_state_t *state = (const tsqr_state_t *) vstate;
  return state->R;
}

static const gsl_vector *
tsqr_QTb(const void * vstate)
{
  const tsqr_state_t *state = (const tsqr_state_t *) vstate;
  return state->QTb;
}

static int
tsqr_rcond(double * rcond, void * vstate)
{
  tsqr_state_t *state = (tsqr_state_t *) vstate;
  return gsl_linalg_tri_rcond(CblasUpper, state->R, rcond, state->work3);
}

/*
tsqr_svd()
  Compute the SVD of the upper triangular
R factor. This allows us to compute the upper/lower
bounds on the regularization parameter and compute
the matrix reciprocal condition number.

Inputs: state - workspace

Return: success/error
*/

static int
tsqr_svd(tsqr_state_t * state)
{
  int status;

  status = gsl_multifit_linear_svd(state->R, state->multifit_workspace_p);
  if (status)
    {
      GSL_ERROR("error computing SVD of R", status);
    }

  state->svd = 1;

  return GSL_SUCCESS;
}

static const gsl_multilarge_linear_type tsqr_type =
{
  "tsqr",
  tsqr_alloc,
  tsqr_reset,
  tsqr_accumulate,
  tsqr_solve,
  tsqr_rcond,
  tsqr_lcurve,
  tsqr_R,
  tsqr_QTb,
  tsqr_free
};

const gsl_multilarge_linear_type * gsl_multilarge_linear_tsqr = &tsqr_type;
