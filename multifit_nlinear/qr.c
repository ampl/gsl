/* multifit_nlinear/qr.c
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

/*
 * This module handles the solution of the linear least squares
 * system:
 *
 * [     J~      ] p~ = - [ f ]
 * [ sqrt(mu)*D~ ]        [ 0 ]
 *
 * using a QR approach. Quantities are scaled according to:
 *
 * J~ = J S
 * D~ = D S
 * p~ = S^{-1} p
 *
 * where S is a diagonal matrix and S_jj = || J_j || and J_j is column
 * j of the Jacobian. This balancing transformation seems to be more
 * numerically stable for some Jacobians.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

#include "common.c"
#include "qrsolv.c"

typedef struct
{
  size_t p;
  gsl_matrix *QR;            /* QR factorization of J */
  gsl_vector *tau_Q;         /* Householder scalars for Q */
  gsl_matrix *T;             /* workspace matrix for qrsolv, p-by-p */
  gsl_permutation *perm;     /* permutation matrix */
  size_t rank;               /* rank of J */
  gsl_vector *residual;      /* residual of LS problem [ J; sqrt(mu) D ] p = - [ f; 0 ] */
  gsl_vector *qtf;           /* Q^T f */
  gsl_vector *workn;         /* workspace, length n */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *work3p;        /* workspace, length 3*p */
  double mu;                 /* LM parameter */
} qr_state_t;

static int qr_init(const void * vtrust_state, void * vstate);
static int qr_presolve(const double mu, const void * vtrust_state, void * vstate);
static int qr_solve(const gsl_vector * f, gsl_vector *x,
                    const void * vtrust_state, void *vstate);
static int qr_rcond(double * rcond, void * vstate);

static void *
qr_alloc (const size_t n, const size_t p)
{
  qr_state_t *state;

  (void)n;
  
  state = calloc(1, sizeof(qr_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate qr state", GSL_ENOMEM);
    }

  state->QR = gsl_matrix_alloc(n, p);
  if (state->QR == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for QR", GSL_ENOMEM);
    }

  state->tau_Q = gsl_vector_alloc(p);
  if (state->tau_Q == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for tau_Q",
                      GSL_ENOMEM);
    }

  state->T = gsl_matrix_alloc(p, p);
  if (state->T == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for T", GSL_ENOMEM);
    }

  state->qtf = gsl_vector_alloc(n);
  if (state->qtf == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for qtf",
                      GSL_ENOMEM);
    }

  state->residual = gsl_vector_alloc(n);
  if (state->residual == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for residual",
                      GSL_ENOMEM);
    }

  state->perm = gsl_permutation_calloc(p);
  if (state->perm == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for perm",
                      GSL_ENOMEM);
    }

  state->workn = gsl_vector_alloc(n);
  if (state->workn == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workn",
                      GSL_ENOMEM);
    }

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workp",
                      GSL_ENOMEM);
    }

  state->work3p = gsl_vector_alloc(3 * p);
  if (state->work3p == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for work3p",
                      GSL_ENOMEM);
    }

  state->p = p;
  state->mu = 0.0;
  state->rank = 0;

  return state;
}

static void
qr_free(void *vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;

  if (state->QR)
    gsl_matrix_free(state->QR);

  if (state->tau_Q)
    gsl_vector_free(state->tau_Q);

  if (state->T)
    gsl_matrix_free(state->T);

  if (state->qtf)
    gsl_vector_free(state->qtf);

  if (state->residual)
    gsl_vector_free(state->residual);

  if (state->perm)
    gsl_permutation_free(state->perm);

  if (state->workn)
    gsl_vector_free(state->workn);

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->work3p)
    gsl_vector_free(state->work3p);

  free(state);
}

/* compute J = Q R PT */
static int
qr_init(const void * vtrust_state, void * vstate)
{
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  qr_state_t *state = (qr_state_t *) vstate;
  int signum;

  /* perform QR decomposition of J */
  gsl_matrix_memcpy(state->QR, trust_state->J);
  gsl_linalg_QRPT_decomp(state->QR, state->tau_Q, state->perm,
                         &signum, state->workp);

  return GSL_SUCCESS;
}

static int
qr_presolve(const double mu, const void * vtrust_state, void * vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;

  state->mu = mu;

  (void) vtrust_state;

  return GSL_SUCCESS;
}

static int
qr_solve(const gsl_vector * f, gsl_vector *x,
         const void * vtrust_state, void *vstate)
{
  qr_state_t *state = (qr_state_t *) vstate;
  int status;

  if (state->mu == 0.0)
    {
      /*
       * compute Gauss-Newton direction by solving
       * J x = f
       * with an attempt to identify rank deficiency in J
       */
      size_t rank = gsl_linalg_QRPT_rank(state->QR, -1.0);
      status = gsl_linalg_QRPT_lssolve2(state->QR, state->tau_Q, state->perm,
                                        f, rank, x, state->residual);
    }
  else
    {
      /*
       * solve:
       *
       * [     J      ] x = [ f ]
       * [ sqrt(mu) D ]     [ 0 ]
       *
       * using QRPT factorization of J
       */

      const gsl_multifit_nlinear_trust_state *trust_state =
        (const gsl_multifit_nlinear_trust_state *) vtrust_state;
      double sqrt_mu = sqrt(state->mu);

      /* compute qtf = Q^T f */
      gsl_vector_memcpy(state->qtf, f);
      gsl_linalg_QR_QTvec(state->QR, state->tau_Q, state->qtf);

      status = qrsolv(state->QR, state->perm, sqrt_mu, trust_state->diag,
                      state->qtf, state->T, x, state->workn);
    }

  /* reverse step to go downhill */
  gsl_vector_scale(x, -1.0);

  return status;
}

static int
qr_rcond(double * rcond, void * vstate)
{
  int status;
  qr_state_t *state = (qr_state_t *) vstate;

  status = gsl_linalg_QRPT_rcond(state->QR, rcond, state->work3p);

  return status;
}

static const gsl_multifit_nlinear_solver qr_type =
{
  "qr",
  qr_alloc,
  qr_init,
  qr_presolve,
  qr_solve,
  qr_rcond,
  qr_free
};

const gsl_multifit_nlinear_solver *gsl_multifit_nlinear_solver_qr = &qr_type;
