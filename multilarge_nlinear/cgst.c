/* multilarge_nlinear/cgst.c
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

/*
 * This module contains an implementation of the Steihaug-Toint
 * conjugate gradient algorithm for nonlinear optimization problems.
 * This implementation closely follows the following works:
 *
 * [1] T. Steihaug, The conjugate gradient method and trust regions
 *     in large scale optimization, SIAM J. Num. Anal., 20(3) 1983.
 *
 * In the below algorithm, the Jacobian and gradient are scaled
 * according to:
 *
 * J~ = J D^{-1}
 * g~ = D^{-1}
 *
 * prior to any calculations which results in better numerical
 * stability when solving for the Gauss-Newton step. The resulting
 * step vector is then backtransformed as:
 *
 * dx = D^{-1} dx~
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *z;             /* Gauss-Newton step, size p */
  gsl_vector *r;             /* steepest descent step, size p */
  gsl_vector *d;             /* steepest descent step, size p */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */
  double norm_g;             /* || g~ || */

  double cgtol;              /* tolerance for CG solution */
  size_t cgmaxit;            /* maximum CG iterations */
} cgst_state_t;

#include "common.c"

static void * cgst_alloc (const void * params, const size_t n, const size_t p);
static void cgst_free(void *vstate);
static int cgst_init(const void *vtrust_state, void *vstate);
static int cgst_preloop(const void * vtrust_state, void * vstate);
static int cgst_step(const void * vtrust_state, const double delta,
                     gsl_vector * dx, void * vstate);
static int cgst_preduction(const void * vtrust_state, const gsl_vector * dx,
                           double * pred, void * vstate);
static double cgst_calc_tau(const gsl_vector * p, const gsl_vector * d,
                            const double delta);

static void *
cgst_alloc (const void * params, const size_t n, const size_t p)
{
  const gsl_multilarge_nlinear_parameters *par = (const gsl_multilarge_nlinear_parameters *) params;
  cgst_state_t *state;
  
  state = calloc(1, sizeof(cgst_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate st state", GSL_ENOMEM);
    }

  state->z = gsl_vector_alloc(p);
  if (state->z == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for z", GSL_ENOMEM);
    }

  state->r = gsl_vector_alloc(p);
  if (state->r == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for r", GSL_ENOMEM);
    }

  state->d = gsl_vector_alloc(p);
  if (state->d == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for d", GSL_ENOMEM);
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

  state->n = n;
  state->p = p;

  state->cgmaxit = par->max_iter;
  if (state->cgmaxit == 0)
    state->cgmaxit = n;

  state->cgtol = par->tol;

  return state;
}

static void
cgst_free(void *vstate)
{
  cgst_state_t *state = (cgst_state_t *) vstate;

  if (state->z)
    gsl_vector_free(state->z);

  if (state->r)
    gsl_vector_free(state->r);

  if (state->d)
    gsl_vector_free(state->d);

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->workn)
    gsl_vector_free(state->workn);

  free(state);
}

/*
cgst_init()
  Initialize cgst solver

Inputs: vtrust_state - trust state
        vstate       - workspace

Return: success/error
*/

static int
cgst_init(const void *vtrust_state, void *vstate)
{
  /* nothing to do */

  (void)vtrust_state;
  (void)vstate;

  return GSL_SUCCESS;
}

static int
cgst_preloop(const void * vtrust_state, void * vstate)
{
  /* nothing to do */

  (void)vtrust_state;
  (void)vstate;

  return GSL_SUCCESS;
}

/*
cgst_step()
  Calculate a new step vector

Return:
GSL_SUCCESS if CG solution found
GSL_EMAXITER if no solution found
*/

static int
cgst_step(const void * vtrust_state, const double delta,
          gsl_vector * dx, void * vstate)
{
  int status;
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  cgst_state_t *state = (cgst_state_t *) vstate;
  const gsl_vector * x = trust_state->x;
  const gsl_vector * f = trust_state->f;
  const gsl_vector * swts = trust_state->sqrt_wts;
  const gsl_vector * diag = trust_state->diag;
  const gsl_multilarge_nlinear_parameters * params = trust_state->params;
  gsl_multilarge_nlinear_fdf * fdf = trust_state->fdf;
  double alpha, beta, u;
  double norm_Jd;   /* || J D^{-1} d_i || */
  double norm_r;    /* || r_i || */
  double norm_rp1;  /* || r_{i+1} || */
  size_t i;

  /* Step 1 of [1], section 2; scale gradient as
   *
   * g~ = D^{-1} g
   *
   * for better numerical stability
   */

  for (i = 0; i < state->p; ++i)
    {
      double gi = gsl_vector_get(trust_state->g, i);
      double di = gsl_vector_get(trust_state->diag, i);

      gsl_vector_set(state->z, i, 0.0);
      gsl_vector_set(state->r, i, -gi / di);
      gsl_vector_set(state->d, i, -gi / di);
      gsl_vector_set(state->workp, i, gi / di);
    }

  /* compute || g~ || */
  state->norm_g = gsl_blas_dnrm2(state->workp);

  for (i = 0; i < state->cgmaxit; ++i)
    {
      /* workp := D^{-1} d_i */
      gsl_vector_memcpy(state->workp, state->d);
      gsl_vector_div(state->workp, trust_state->diag);

      /* workn := J D^{-1} d_i */
      status = gsl_multilarge_nlinear_eval_df(CblasNoTrans, x, f, state->workp,
                                              swts, params->h_df, params->fdtype,
                                              fdf, state->workn, NULL, NULL);
      if (status)
        return status;

      /* compute || J D^{-1} d_i || */
      norm_Jd = gsl_blas_dnrm2(state->workn);

      /* Step 2 of [1], section 2 */
      if (norm_Jd == 0.0)
        {
          double tau = cgst_calc_tau(state->z, state->d, delta);

          /* dx = z_i + tau*d_i */
          scaled_addition(1.0, state->z, tau, state->d, dx);
          gsl_vector_div(dx, diag);

          return GSL_SUCCESS;
        }

      /* Step 3 of [1], section 2 */

      norm_r = gsl_blas_dnrm2(state->r);
      u = norm_r / norm_Jd;
      alpha = u * u;

      /* workp <= z_{i+1} = z_i + alpha_i*d_i */
      scaled_addition(1.0, state->z, alpha, state->d, state->workp);

      u = gsl_blas_dnrm2(state->workp);
      if (u >= delta)
        {
          double tau = cgst_calc_tau(state->z, state->d, delta);

          /* dx = z_i + tau*d_i */
          scaled_addition(1.0, state->z, tau, state->d, dx);
          gsl_vector_div(dx, diag);

          return GSL_SUCCESS;
        }

      /* store z_{i+1} */
      gsl_vector_memcpy(state->z, state->workp);

      /* Step 4 of [1], section 2 */

      /* compute: workp := alpha B d_i = alpha D^{-1} J^T J D^{-1} d_i,
       * where J D^{-1} d_i is already stored in workn */
      status = gsl_multilarge_nlinear_eval_df(CblasTrans, x, f, state->workn,
                                              swts, params->h_df, params->fdtype,
                                              fdf, state->workp, NULL, NULL);
      if (status)
        return status;

      gsl_vector_div(state->workp, trust_state->diag);
      gsl_vector_scale(state->workp, alpha);

      /* r_{i+1} = r_i - alpha*B*d_i */
      gsl_vector_sub(state->r, state->workp);
      norm_rp1 = gsl_blas_dnrm2(state->r);

      u = norm_rp1 / state->norm_g;
      if (u < state->cgtol)
        {
          gsl_vector_memcpy(dx, state->z);
          gsl_vector_div(dx, diag);
          return GSL_SUCCESS;
        }

      /* Step 5 of [1], section 2 */

      /* compute u = ||r_{i+1}|| / || r_i|| */
      u = norm_rp1 / norm_r;
      beta = u * u;

      /* compute: d_{i+1} = rt_{i+1} + beta*d_i */
      scaled_addition(1.0, state->r, beta, state->d, state->d);
    }

  /* failed to converge, return current estimate */
  gsl_vector_memcpy(dx, state->z);
  gsl_vector_div(dx, diag);

  return GSL_EMAXITER;
}

static int
cgst_preduction(const void * vtrust_state, const gsl_vector * dx,
                double * pred, void * vstate)
{
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  cgst_state_t *state = (cgst_state_t *) vstate;

  *pred = quadratic_preduction(trust_state, dx, state->workn);

  return GSL_SUCCESS;
}

/*
cgst_calc_tau()
  Compute tau > 0 such that:

|| p + tau*d || = delta
*/

static double
cgst_calc_tau(const gsl_vector * p, const gsl_vector * d,
              const double delta)
{
  double norm_p, norm_d, u;
  double t1, t2, tau;

  norm_p = gsl_blas_dnrm2(p);
  norm_d = gsl_blas_dnrm2(d);

  /* compute (p, d) */
  gsl_blas_ddot(p, d, &u);

  t1 = u / (norm_d * norm_d);
  t2 = t1*u + (delta + norm_p) * (delta - norm_p);
  tau = -t1 + sqrt(t2) / norm_d;

  return tau;
}

static const gsl_multilarge_nlinear_trs cgst_type =
{
  "steihaug-toint",
  cgst_alloc,
  cgst_init,
  cgst_preloop,
  cgst_step,
  cgst_preduction,
  cgst_free
};

const gsl_multilarge_nlinear_trs *gsl_multilarge_nlinear_trs_cgst = &cgst_type;
