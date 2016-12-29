/* multilarge_nlinear/dogleg.c
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
 * This module contains an implementation of the Powell dogleg
 * algorithm for nonlinear optimization problems. This implementation
 * closely follows the following works:
 *
 * [1] H. B. Nielsen, K. Madsen, Introduction to Optimization and
 *     Data Fitting, Informatics and Mathematical Modeling,
 *     Technical University of Denmark (DTU), 2010.
 *
 * [2] J. E. Dennis and H. H. W. Mei, Two new unconstrained optimization
 *     algorithms which use function and gradient values, J. Opt. Theory and
 *     Appl., 28(4), 1979.
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *dx_gn;         /* Gauss-Newton step, size p */
  gsl_vector *dx_sd;         /* steepest descent step, size p */
  double norm_Dgn;           /* || D dx_gn || */
  double norm_Dsd;           /* || D dx_sd || */
  double norm_Dinvg;         /* || D^{-1} g || */
  double norm_JDinv2g;       /* || J D^{-2} g || */
  gsl_vector *workp1;        /* workspace, length p */
  gsl_vector *workp2;        /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */

  /* tunable parameters */
  gsl_multilarge_nlinear_parameters params;
} dogleg_state_t;

#include "common.c"

static void * dogleg_alloc (const void * params, const size_t n, const size_t p);
static void dogleg_free(void *vstate);
static int dogleg_init(const void *vtrust_state, void *vstate);
static int dogleg_preloop(const void * vtrust_state, void * vstate);
static int dogleg_step(const void * vtrust_state, const double delta,
                       gsl_vector * dx, void * vstate);
static int dogleg_double_step(const void * vtrust_state, const double delta,
                              gsl_vector * dx, void * vstate);
static int dogleg_preduction(const void * vtrust_state, const gsl_vector * dx,
                             double * pred, void * vstate);
static int dogleg_calc_gn(const gsl_multilarge_nlinear_trust_state * trust_state, gsl_vector * dx);
static double dogleg_beta(const double t, const double delta,
                          const gsl_vector * diag, dogleg_state_t * state);

static void *
dogleg_alloc (const void * params, const size_t n, const size_t p)
{
  dogleg_state_t *state;
  
  state = calloc(1, sizeof(dogleg_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate dogleg state", GSL_ENOMEM);
    }

  state->dx_gn = gsl_vector_alloc(p);
  if (state->dx_gn == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for dx_gn", GSL_ENOMEM);
    }

  state->dx_sd = gsl_vector_alloc(p);
  if (state->dx_sd == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for dx_sd", GSL_ENOMEM);
    }

  state->workp1 = gsl_vector_alloc(p);
  if (state->workp1 == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workp1", GSL_ENOMEM);
    }

  state->workp2 = gsl_vector_alloc(p);
  if (state->workp2 == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workp2", GSL_ENOMEM);
    }

  state->workn = gsl_vector_alloc(n);
  if (state->workn == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workn", GSL_ENOMEM);
    }

  state->n = n;
  state->p = p;
  state->params = *(const gsl_multilarge_nlinear_parameters *) params;

  return state;
}

static void
dogleg_free(void *vstate)
{
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  if (state->dx_gn)
    gsl_vector_free(state->dx_gn);

  if (state->dx_sd)
    gsl_vector_free(state->dx_sd);

  if (state->workp1)
    gsl_vector_free(state->workp1);

  if (state->workp2)
    gsl_vector_free(state->workp2);

  if (state->workn)
    gsl_vector_free(state->workn);

  free(state);
}

/*
dogleg_init()
  Initialize dogleg solver

Inputs: vtrust_state - trust state
        vstate       - workspace

Return: success/error
*/

static int
dogleg_init(const void *vtrust_state, void *vstate)
{
  (void)vtrust_state;
  (void)vstate;

  return GSL_SUCCESS;
}

/*
dogleg_preloop()
  Initialize dogleg method prior to iteration loop.
This involves computing the steepest descent step. The
Gauss-Newton step is computed later in the _step() functions
if required.

Notes: on output,
1) state->dx_sd contains steepest descent step
2) state->norm_Dinvg contains || D^{-1} g ||
3) state->norm_JDinv2g contains || J D^{-2} g ||
*/

static int
dogleg_preloop(const void * vtrust_state, void * vstate)
{
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;
  double u;
  double alpha; /* ||g||^2 / ||Jg||^2 */

  /* calculate the steepest descent step */

  /* compute workp1 = D^{-1} g and its norm */
  gsl_vector_memcpy(state->workp1, trust_state->g);
  gsl_vector_div(state->workp1, trust_state->diag);
  state->norm_Dinvg = gsl_blas_dnrm2(state->workp1);

  /* compute workp1 = D^{-2} g */
  gsl_vector_div(state->workp1, trust_state->diag);

  /* compute workp2 = J^T J D^{-2} g */
  gsl_blas_dsymv(CblasLower, 1.0, trust_state->JTJ, state->workp1, 0.0, state->workp2);

  /* compute norm_JDinv2g = || J D^{-2} g || */
  gsl_blas_ddot(state->workp1, state->workp2, &u);
  state->norm_JDinv2g = sqrt(u);

  u = state->norm_Dinvg / state->norm_JDinv2g;
  alpha = u * u;

  /* dx_sd = -alpha D^{-2} g */
  gsl_vector_memcpy(state->dx_sd, state->workp1);
  gsl_vector_scale(state->dx_sd, -alpha);

  state->norm_Dsd = scaled_enorm(trust_state->diag, state->dx_sd);
  state->norm_Dgn = -1.0; /* computed later if needed */

  return GSL_SUCCESS;
}

/*
dogleg_step()
  Calculate a new step vector
*/

static int
dogleg_step(const void * vtrust_state, const double delta,
            gsl_vector * dx, void * vstate)
{
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  if (state->norm_Dsd >= delta)
    {
      /* steepest descent step is outside trust region;
       * truncate steepest descent step to trust region boundary */
      gsl_vector_memcpy(dx, state->dx_sd);
      gsl_vector_scale(dx, delta / state->norm_Dsd);
    }
  else
    {
      /* compute Gauss-Newton step if needed */
      if (state->norm_Dgn < 0.0)
        {
          int status = dogleg_calc_gn(trust_state, state->dx_gn);

          if (status)
            return status;

          /* compute || D dx_gn || */
          state->norm_Dgn = scaled_enorm(trust_state->diag, state->dx_gn);
        }

      if (state->norm_Dgn <= delta)
        {
          /* Gauss-Newton step is inside trust region, use it as final step
           * since it is the global minimizer of the quadratic model function */
          gsl_vector_memcpy(dx, state->dx_gn);
        }
      else
        {
          /* Gauss-Newton step is outside trust region, but steepest
           * descent is inside; use dogleg step */

          double beta = dogleg_beta(1.0, delta, trust_state->diag, state);

          /* compute: workp1 = dx_gn - dx_sd */
          scaled_addition(1.0, state->dx_gn, -1.0, state->dx_sd, state->workp1);

          /* dx = dx_sd + beta*(dx_gn - dx_sd) */
          scaled_addition(beta, state->workp1, 1.0, state->dx_sd, dx);
        }
    }

  return GSL_SUCCESS;
}

/*
dogleg_double_step()
  Calculate a new step with double dogleg method. Based on
section 3 of [2]
*/

static int
dogleg_double_step(const void * vtrust_state, const double delta,
                   gsl_vector * dx, void * vstate)
{
  const double alpha_fac = 0.8; /* recommended value from Dennis and Mei */
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  if (state->norm_Dsd >= delta)
    {
      /* steepest descent step is outside trust region;
       * truncate steepest descent step to trust region boundary */
      gsl_vector_memcpy(dx, state->dx_sd);
      gsl_vector_scale(dx, delta / state->norm_Dsd);
    }
  else
    {
      /* compute Gauss-Newton step if needed */
      if (state->norm_Dgn < 0.0)
        {
          int status = dogleg_calc_gn(trust_state, state->dx_gn);

          if (status)
            return status;

          /* compute || D dx_gn || */
          state->norm_Dgn = scaled_enorm(trust_state->diag, state->dx_gn);
        }

      if (state->norm_Dgn <= delta)
        {
          /* Gauss-Newton step is inside trust region, use it as final step
           * since it is the global minimizer of the quadratic model function */
          gsl_vector_memcpy(dx, state->dx_gn);
        }
      else
        {
          double t, u, v, c;

          /* compute: u = ||D^{-1} g||^2 / ||J D^{-2} g||^2 */
          v = state->norm_Dinvg / state->norm_JDinv2g;
          u = v * v;

          /* compute: v = g^T dx_gn */
          gsl_blas_ddot(trust_state->g, state->dx_gn, &v);

          /* compute: c = ||D^{-1} g||^4 / (||J D^{-2} g||^2 * |g^T dx_gn|) */
          c = u * (state->norm_Dinvg / fabs(v)) * state->norm_Dinvg;

          /* compute: t = 1 - alpha_fac*(1-c) */
          t = 1.0 - alpha_fac*(1.0 - c);

          if (t * state->norm_Dgn <= delta)
            {
              /* set dx = (delta / ||D dx_gn||) dx_gn */
              gsl_vector_memcpy(dx, state->dx_gn);
              gsl_vector_scale(dx, delta / state->norm_Dgn);
            }
          else
            {
              /* Cauchy point is inside, Gauss-Newton is outside trust region;
               * use double dogleg step */

              double beta = dogleg_beta(t, delta, trust_state->diag, state);

              /* compute: workp1 = t*dx_gn - dx_sd */
              scaled_addition(t, state->dx_gn, -1.0, state->dx_sd, state->workp1);

              /* dx = dx_sd + beta*(t*dx_gn - dx_sd) */
              scaled_addition(beta, state->workp1, 1.0, state->dx_sd, dx);
            }
        }
    }

  return GSL_SUCCESS;
}

static int
dogleg_preduction(const void * vtrust_state, const gsl_vector * dx,
                  double * pred, void * vstate)
{
  const gsl_multilarge_nlinear_trust_state *trust_state =
    (const gsl_multilarge_nlinear_trust_state *) vtrust_state;
  dogleg_state_t *state = (dogleg_state_t *) vstate;

  *pred = quadratic_preduction(trust_state, dx, state->workn);

  return GSL_SUCCESS;
}

/*
dogleg_calc_gn()
  Calculate Gauss-Newton step by solving

J^T J dx_gn = -J^T f

Inputs: trust_state - trust state variables
        dx          - (output) Gauss-Newton step vector

Return: success/error
*/

static int
dogleg_calc_gn(const gsl_multilarge_nlinear_trust_state * trust_state, gsl_vector * dx)
{
  int status;
  const gsl_multilarge_nlinear_parameters *params = trust_state->params;

  /* initialize linear least squares solver */
  status = (params->solver->init)(trust_state, trust_state->solver_state);
  if (status)
    return status;

  /* prepare the linear solver to compute Gauss-Newton step */
  status = (params->solver->presolve)(0.0, trust_state, trust_state->solver_state);
  if (status)
    return status;

  /* solve: J dx_gn = -f for Gauss-Newton step */
  status = (params->solver->solve)(trust_state->g,
                                   dx,
                                   trust_state,
                                   trust_state->solver_state);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
dogleg_beta()
  This function finds beta in [0,1] such that the step

dx = dx_sd + beta*(t*dx_gn - dx_sd)

has norm

||D dx|| = delta

beta is the positive root of the quadratic:

a beta^2 + b beta + c = 0

with

a = ||D(t*dx_gn - dx_sd)||^2
b = 2 dx_sd^T D^T D (t*dx_gn - dx_sd)
c = ||D dx_sd||^2 - delta^2 

Inputs: t     - amount of Gauss-Newton step to use for dogleg
                (= 1 for classical dogleg, <= 1 for double dogleg)
        delta - trust region radius
        diag  - diag(D) scaling matrix
        state - workspace
*/

static double
dogleg_beta(const double t, const double delta,
            const gsl_vector * diag, dogleg_state_t * state)
{
  double beta;
  double a, b, c;

  /* compute: workp1 = t*dx_gn - dx_sd */
  scaled_addition(t, state->dx_gn, -1.0, state->dx_sd, state->workp1);

  /* a = || D (t*dx_gn - dx_sd) ||^2 */
  a = scaled_enorm(diag, state->workp1);
  a *= a;

  /* workp1 = D^T D (t*dx_gn - dx_sd) */
  gsl_vector_mul(state->workp1, diag);
  gsl_vector_mul(state->workp1, diag);

  /* b = 2 dx_sd^T D^T D (t*dx_gn - dx-sd) */
  gsl_blas_ddot(state->dx_sd, state->workp1, &b);
  b *= 2.0;

  /* c = || D dx_sd ||^2 - delta^2 = (||D dx_sd|| + delta) (||D dx_sd|| - delta) */
  c = (state->norm_Dsd + delta) * (state->norm_Dsd - delta);

  if (b > 0.0)
    {
      beta = (-2.0 * c) / (b + sqrt(b*b - 4.0*a*c));
    }
  else
    {
      beta = (-b + sqrt(b*b - 4.0*a*c)) / (2.0 * a);
    }

  return beta;
}

static const gsl_multilarge_nlinear_trs dogleg_type =
{
  "dogleg",
  dogleg_alloc,
  dogleg_init,
  dogleg_preloop,
  dogleg_step,
  dogleg_preduction,
  dogleg_free
};

const gsl_multilarge_nlinear_trs *gsl_multilarge_nlinear_trs_dogleg = &dogleg_type;

static const gsl_multilarge_nlinear_trs ddogleg_type =
{
  "double-dogleg",
  dogleg_alloc,
  dogleg_init,
  dogleg_preloop,
  dogleg_double_step,
  dogleg_preduction,
  dogleg_free
};

const gsl_multilarge_nlinear_trs *gsl_multilarge_nlinear_trs_ddogleg = &ddogleg_type;
