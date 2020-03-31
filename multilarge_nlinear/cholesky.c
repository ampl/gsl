/* multilarge_nlinear/cholesky.c
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
 * This module calculates the solution of the normal equations least squares
 * system:
 *
 * [ J^T J + mu D^T D ] p = -J^T f
 *
 * using the Cholesky decomposition.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

#include "common.c"

typedef struct
{
  gsl_matrix *JTJ;           /* J^T J */
  gsl_matrix *work_JTJ;      /* copy of J^T J */
  gsl_vector *rhs;           /* -J^T f, size p */
  gsl_vector *work3p;        /* workspace, size 3*p */
  gsl_vector *workn;         /* workspace, size n */
  double mu;                 /* current regularization parameter */
} cholesky_state_t;

static void *cholesky_alloc (const size_t n, const size_t p);
static int cholesky_init(const void * vtrust_state, void * vstate);
static int cholesky_presolve(const double mu, const void * vtrust_state, void * vstate);
static int cholesky_solve(const gsl_vector * g, gsl_vector *x,
                          const  void * vtrust_state, void *vstate);
static int cholesky_rcond(double * rcond, const gsl_matrix * JTJ, void * vstate);
static int cholesky_covar(const gsl_matrix * JTJ, gsl_matrix * covar, void * vstate);
static int cholesky_solve_rhs(const gsl_vector * b, gsl_vector *x, cholesky_state_t *state);
static int cholesky_regularize(const double mu, const gsl_vector * diag, gsl_matrix * A,
                               cholesky_state_t * state);

static void *
cholesky_alloc (const size_t n, const size_t p)
{
  cholesky_state_t *state;

  state = calloc(1, sizeof(cholesky_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate cholesky state", GSL_ENOMEM);
    }

  state->JTJ = gsl_matrix_alloc(p, p);
  if (state->JTJ == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for JTJ", GSL_ENOMEM);
    }

  state->work_JTJ = gsl_matrix_alloc(p, p);
  if (state->work_JTJ == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for JTJ workspace",
                      GSL_ENOMEM);
    }

  state->rhs = gsl_vector_alloc(p);
  if (state->rhs == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for rhs", GSL_ENOMEM);
    }

  state->work3p = gsl_vector_alloc(3 * p);
  if (state->work3p == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for work3p", GSL_ENOMEM);
    }

  state->workn = gsl_vector_alloc(n);
  if (state->workn == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workn", GSL_ENOMEM);
    }

  state->mu = -1.0;

  return state;
}

static void
cholesky_free(void *vstate)
{
  cholesky_state_t *state = (cholesky_state_t *) vstate;

  if (state->JTJ)
    gsl_matrix_free(state->JTJ);

  if (state->work_JTJ)
    gsl_matrix_free(state->work_JTJ);

  if (state->rhs)
    gsl_vector_free(state->rhs);

  if (state->work3p)
    gsl_vector_free(state->work3p);

  if (state->workn)
    gsl_vector_free(state->workn);

  free(state);
}

static int
cholesky_init(const void * vtrust_state, void * vstate)
{
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  cholesky_state_t *state = (cholesky_state_t *) vstate;

  /* store J^T J normal equations matrix */
  gsl_matrix_tricpy(CblasLower, CblasNonUnit, state->JTJ, trust_state->JTJ);

  return GSL_SUCCESS;
}

/*
cholesky_presolve()
  Compute the modified Cholesky decomposition of J^T J + mu D^T D.
Modified Cholesky is used in case mu = 0 and there are rounding
errors in forming J^T J which could lead to an indefinite matrix.

Inputs: mu     - LM parameter
        vstate - workspace

Notes:
1) On output, state->work_JTJ contains the Cholesky decomposition of
J^T J + mu D^T D
*/

static int
cholesky_presolve(const double mu, const void * vtrust_state, void * vstate)
{
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  cholesky_state_t *state = (cholesky_state_t *) vstate;
  gsl_matrix *JTJ = state->work_JTJ;
  const gsl_vector *diag = trust_state->diag;
  int status;

  /* copy lower triangle of A to workspace */
  gsl_matrix_tricpy(CblasLower, CblasNonUnit, JTJ, state->JTJ);

  /* augment normal equations: A -> A + mu D^T D */
  status = cholesky_regularize(mu, diag, JTJ, state);
  if (status)
    return status;

  /* compute Cholesky decomposition */
  status = gsl_linalg_cholesky_decomp1(JTJ);
  if (status)
    return status;

  state->mu = mu;

  return GSL_SUCCESS;
}

/*
cholesky_solve()
  Compute (J^T J + mu D^T D) x = -g

where g = J^T f

Inputs: g      - right hand side vector g, size p
        x      - (output) solution vector
        vstate - cholesky workspace
*/

static int
cholesky_solve(const gsl_vector * g, gsl_vector *x,
               const void * vtrust_state, void *vstate)
{
  cholesky_state_t *state = (cholesky_state_t *) vstate;
  int status;

  status = cholesky_solve_rhs(g, x, state);
  if (status)
    return status;

  /* reverse direction to go downhill */
  gsl_vector_scale(x, -1.0);

  (void) vtrust_state;

  return GSL_SUCCESS;
}

static int
cholesky_rcond(double * rcond, const gsl_matrix * JTJ, void * vstate)
{
  int status;
  cholesky_state_t *state = (cholesky_state_t *) vstate;
  gsl_error_handler_t * err_handler;
  double rcond_JTJ;

  /* its possible the current Cholesky decomposition is from the previous
   * iteration so do a new one to be sure we use the right Jacobian */

  /* copy lower triangle of JTJ to workspace */
  gsl_matrix_tricpy(CblasLower, CblasNonUnit, state->work_JTJ, JTJ);

  /* compute Cholesky decomposition, turning off error handler */
  err_handler = gsl_set_error_handler_off();
  status = gsl_linalg_cholesky_decomp1(state->work_JTJ);
  gsl_set_error_handler(err_handler);

  if (status)
    {
      /* matrix is singular */
      *rcond = 0.0;
      return GSL_SUCCESS;
    }

  status = gsl_linalg_cholesky_rcond(state->work_JTJ, &rcond_JTJ, state->work3p);
  if (status == GSL_SUCCESS)
    *rcond = sqrt(rcond_JTJ);

  return status;
}

static int
cholesky_covar(const gsl_matrix * JTJ, gsl_matrix * covar, void * vstate)
{
  int status;
  gsl_error_handler_t * err_handler;

  (void) vstate;

  /* its possible the current Cholesky decomposition is from the previous
   * iteration so do a new one to be sure we use the right Jacobian */

  /* copy lower triangle of JTJ to workspace */
  gsl_matrix_tricpy(CblasLower, CblasNonUnit, covar, JTJ);

  /* compute Cholesky decomposition, turning off error handler */
  err_handler = gsl_set_error_handler_off();
  status = gsl_linalg_cholesky_decomp1(covar);
  gsl_set_error_handler(err_handler);

  if (status)
    return status;

  status = gsl_linalg_cholesky_invert(covar);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/* solve: (J^T J + mu D^T D) x = b */
static int
cholesky_solve_rhs(const gsl_vector * b, gsl_vector *x, cholesky_state_t *state)
{
  int status;
  gsl_matrix *JTJ = state->work_JTJ;

  status = gsl_linalg_cholesky_solve(JTJ, b, x);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/* A <- A + mu D^T D */
static int
cholesky_regularize(const double mu, const gsl_vector * diag, gsl_matrix * A,
                    cholesky_state_t * state)
{
  (void) state;

  if (mu != 0.0)
    {
      size_t i;

      for (i = 0; i < diag->size; ++i)
        {
          double di = gsl_vector_get(diag, i);
          double *Aii = gsl_matrix_ptr(A, i, i);
          *Aii += mu * di * di;
        }
    }

  return GSL_SUCCESS;
}

static const gsl_multilarge_nlinear_solver cholesky_type =
{
  "cholesky",
  cholesky_alloc,
  cholesky_init,
  cholesky_presolve,
  cholesky_solve,
  cholesky_rcond,
  cholesky_covar,
  cholesky_free
};

const gsl_multilarge_nlinear_solver *gsl_multilarge_nlinear_solver_cholesky = &cholesky_type;
