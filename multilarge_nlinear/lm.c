/* multilarge_nlinear/lm.c
 * 
 * Copyright (C) 2014, 2015, 2016 Patrick Alken
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
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>

/*
 * This module contains an implementation of the Levenberg-Marquardt
 * algorithm for nonlinear optimization problems. This implementation
 * closely follows the following works:
 *
 * [1] H. B. Nielsen, K. Madsen, Introduction to Optimization and
 *     Data Fitting, Informatics and Mathematical Modeling,
 *     Technical University of Denmark (DTU), 2010.
 *
 * [2] J. J. More, The Levenberg-Marquardt Algorithm: Implementation
 *     and Theory, Lecture Notes in Mathematics, v630, 1978.
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *fvv;           /* D_v^2 f(x), size n */
  gsl_vector *vel;           /* geodesic velocity (standard LM step), size p */
  gsl_vector *acc;           /* geodesic acceleration, size p */
  gsl_vector *JTfvv;         /* J^T fvv, size p */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */

  int accel;                 /* use geodesic acceleration? */

  /* tunable parameters */
  gsl_multilarge_nlinear_parameters params;
} lm_state_t;

#include "common.c"

static void *lm_alloc (const int accel, const void * params, const size_t n, const size_t p);
static void *lm_alloc_noaccel (const void * params, const size_t n, const size_t p);
static void *lm_alloc_accel (const void * params, const size_t n, const size_t p);
static void lm_free(void *vstate);
static int lm_init(const void *vtrust_state, void *vstate);
static int lm_preloop(const void * vtrust_state, void * vstate);
static int lm_step(const void * vtrust_state, const double delta,
                   gsl_vector * dx, void * vstate);
static int lm_preduction(const void * vtrust_state, const gsl_vector * dx,
                         double * pred, void * vstate);

static void *
lm_alloc (const int accel, const void * params, const size_t n, const size_t p)
{
  const gsl_multilarge_nlinear_parameters *mparams = (const gsl_multilarge_nlinear_parameters *) params;
  lm_state_t *state;
  
  state = calloc(1, sizeof(lm_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate lm state", GSL_ENOMEM);
    }

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workp", GSL_ENOMEM);
    }

  state->workn = gsl_vector_alloc(n);
  if (state->workn == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workn", GSL_ENOMEM);
    }

  state->fvv = gsl_vector_alloc(n);
  if (state->fvv == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for fvv", GSL_ENOMEM);
    }

  state->vel = gsl_vector_alloc(p);
  if (state->vel == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for vel", GSL_ENOMEM);
    }

  state->acc = gsl_vector_alloc(p);
  if (state->acc == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for acc", GSL_ENOMEM);
    }

  state->JTfvv = gsl_vector_alloc(p);
  if (state->JTfvv == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for JTfvv", GSL_ENOMEM);
    }

  state->n = n;
  state->p = p;
  state->params = *mparams;
  state->accel = accel;

  return state;
}

static void *
lm_alloc_noaccel (const void * params, const size_t n, const size_t p)
{
  return lm_alloc(0, params, n, p);
}

static void *
lm_alloc_accel (const void * params, const size_t n, const size_t p)
{
  return lm_alloc(1, params, n, p);
}

static void
lm_free(void *vstate)
{
  lm_state_t *state = (lm_state_t *) vstate;

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->workn)
    gsl_vector_free(state->workn);

  if (state->fvv)
    gsl_vector_free(state->fvv);

  if (state->vel)
    gsl_vector_free(state->vel);

  if (state->acc)
    gsl_vector_free(state->acc);

  if (state->JTfvv)
    gsl_vector_free(state->JTfvv);

  free(state);
}

/*
lm_init()
  Initialize LM solver

Inputs: vtrust_state - trust state
        vstate       - workspace

Return: success/error
*/

static int
lm_init(const void *vtrust_state, void *vstate)
{
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  lm_state_t *state = (lm_state_t *) vstate;

  gsl_vector_set_zero(state->vel);
  gsl_vector_set_zero(state->acc);

  *(trust_state->avratio) = 0.0;

  return GSL_SUCCESS;
}

/*
lm_preloop()
  Initialize LM method for new Jacobian matrix
*/

static int
lm_preloop(const void * vtrust_state, void * vstate)
{
  int status;
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  const gsl_multilarge_nlinear_parameters *params = trust_state->params;

  (void)vstate;

  /* initialize linear least squares solver */
  status = (params->solver->init)(trust_state, trust_state->solver_state);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
lm_step()
  Calculate a new step vector by solving the linear
least squares system:

[      J     ] v = - [ f ]
[ sqrt(mu) D ]       [ 0 ]
*/

static int
lm_step(const void * vtrust_state, const double delta,
        gsl_vector * dx, void * vstate)
{
  int status;
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  lm_state_t *state = (lm_state_t *) vstate;
  const gsl_multilarge_nlinear_parameters *params = trust_state->params;
  const double mu = *(trust_state->mu);

  (void)delta;

  /* prepare the linear solver with current LM parameter mu */
  status = (params->solver->presolve)(mu, trust_state, trust_state->solver_state);
  if (status)
    return status;

  /* solve: (J^T J + mu D^T D) v = - J^T f */
  status = (params->solver->solve)(trust_state->g,
                                   state->vel,
                                   trust_state,
                                   trust_state->solver_state);
  if (status)
    return status;

  if (state->accel)
    {
      double anorm, vnorm;

      /* compute geodesic acceleration */
      status = gsl_multilarge_nlinear_eval_fvv(params->h_fvv,
                                               trust_state->x,
                                               state->vel,
                                               trust_state->f,
                                               trust_state->sqrt_wts,
                                               trust_state->fdf,
                                               state->fvv,
                                               state->workp);
      if (status)
        return status;

      /* compute J^T fvv */
      status = gsl_multilarge_nlinear_eval_df(CblasTrans,
                                              trust_state->x,
                                              trust_state->f,
                                              state->fvv,
                                              trust_state->sqrt_wts,
                                              params->h_df,
                                              params->fdtype,
                                              trust_state->fdf,
                                              state->JTfvv,
                                              NULL,
                                              state->workn);
      if (status)
        return status;

      /* solve: (J^T J + mu D^T D) a = - J^T fvv */
      status = (params->solver->solve)(state->JTfvv,
                                       state->acc,
                                       trust_state,
                                       trust_state->solver_state);
      if (status)
        return status;

      anorm = gsl_blas_dnrm2(state->acc);
      vnorm = gsl_blas_dnrm2(state->vel);

      /* store |a| / |v| */
      *(trust_state->avratio) = anorm / vnorm;
    }

  /* compute step dx = v + 1/2 a */
  scaled_addition(1.0, state->vel, 0.5, state->acc, dx);

  return GSL_SUCCESS;
}

/*
lm_preduction()
  Compute predicted reduction using Eq 4.4 of More 1978
*/

static int
lm_preduction(const void * vtrust_state, const gsl_vector * dx,
              double * pred, void * vstate)
{
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  lm_state_t *state = (lm_state_t *) vstate;
  const gsl_vector *diag = trust_state->diag;
  const gsl_vector *p = state->vel;
  const double norm_Dp = scaled_enorm(diag, p);
  const double normf = gsl_blas_dnrm2(trust_state->f);
  const double mu = *(trust_state->mu);
  double u, v;

  (void)dx;

  /* compute workp = J^T J p */
  gsl_blas_dsymv(CblasLower, 1.0, trust_state->JTJ, p, 0.0, state->workp);

  /* compute u = p^T J^T J p = ||J p||^2 */
  gsl_blas_ddot(state->workp, p, &u);

  /* compute u = ||J p||^2 / ||f||^2 */
  u /= normf * normf;

  /* v = ||D p|| / ||f|| */
  v = norm_Dp / normf;

  *pred = u + 2.0 * mu * v * v;

  return GSL_SUCCESS;
}

static const gsl_multilarge_nlinear_trs lm_type =
{
  "levenberg-marquardt",
  lm_alloc_noaccel,
  lm_init,
  lm_preloop,
  lm_step,
  lm_preduction,
  lm_free
};

const gsl_multilarge_nlinear_trs *gsl_multilarge_nlinear_trs_lm = &lm_type;

static const gsl_multilarge_nlinear_trs lmaccel_type =
{
  "levenberg-marquardt+accel",
  lm_alloc_accel,
  lm_init,
  lm_preloop,
  lm_step,
  lm_preduction,
  lm_free
};

const gsl_multilarge_nlinear_trs *gsl_multilarge_nlinear_trs_lmaccel = &lmaccel_type;
