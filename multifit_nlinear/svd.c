/* multifit_nlinear/svd.c
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
 * This module handles the solution of the linear least squares
 * system:
 *
 * [     J      ] dx = - [ f ]
 * [ sqrt(mu)*D ]        [ 0 ]
 *
 * using an SVD approach. The system above is transformed to "standard form"
 * via:
 *
 * J~ = J D^{-1}
 * dx~ = D dx
 *
 * so that
 *
 * [     J~     ] dx~ = - [ f ]
 * [ sqrt(mu)*I ]         [ 0 ]
 *
 * can be solved with a standard SVD method, and then dx is recovered
 * from dx~ via: dx = D^{-1} dx~
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

typedef struct
{
  size_t n;                  /* number of residuals */
  size_t p;                  /* number of parameters */
  gsl_matrix *U;             /* U factor of J, n-by-p */
  gsl_matrix *V;             /* V factor of J, p-by-p */
  gsl_vector *S;             /* singular values, size p */
  gsl_vector *workp;         /* workspace, length p */
  double mu;                 /* LM parameter */
} svd_state_t;

static int svd_init(const void * vtrust_state, void * vstate);
static int svd_presolve(const double mu, const void * vtrust_state, void * vstate);
static int svd_solve(const gsl_vector * f, gsl_vector *x,
                     const void * vtrust_state, void *vstate);
static int svd_rcond(double * rcond, void * vstate);

static void *
svd_alloc (const size_t n, const size_t p)
{
  svd_state_t *state;

  (void)n;
  
  state = calloc(1, sizeof(svd_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate svd state", GSL_ENOMEM);
    }

  state->U = gsl_matrix_alloc(n, p);
  if (state->U == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for U", GSL_ENOMEM);
    }

  state->V = gsl_matrix_alloc(p, p);
  if (state->V == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for V", GSL_ENOMEM);
    }

  state->S = gsl_vector_alloc(p);
  if (state->S == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for S",
                      GSL_ENOMEM);
    }

  state->workp = gsl_vector_alloc(p);
  if (state->workp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for workp",
                      GSL_ENOMEM);
    }

  state->mu = 0.0;
  state->n = n;
  state->p = p;

  return state;
}

static void
svd_free(void *vstate)
{
  svd_state_t *state = (svd_state_t *) vstate;

  if (state->U)
    gsl_matrix_free(state->U);

  if (state->V)
    gsl_matrix_free(state->V);

  if (state->S)
    gsl_vector_free(state->S);

  if (state->workp)
    gsl_vector_free(state->workp);

  free(state);
}

/* compute svd of J */
static int
svd_init(const void * vtrust_state, void * vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  svd_state_t *state = (svd_state_t *) vstate;
  size_t i;

  gsl_matrix_set_zero(state->U);

  /* compute U = J D^{-1} */
  for (i = 0; i < state->p; ++i)
    {
      gsl_vector_const_view Ji = gsl_matrix_const_column(trust_state->J, i);
      gsl_vector_view ui = gsl_matrix_column(state->U, i);
      double di = gsl_vector_get(trust_state->diag, i);

      gsl_blas_daxpy(1.0 / di, &Ji.vector, &ui.vector);
    }

  status = gsl_linalg_SV_decomp(state->U, state->V, state->S, state->workp);

  return status;
}

static int
svd_presolve(const double mu, const void * vtrust_state, void * vstate)
{
  svd_state_t *state = (svd_state_t *) vstate;

  state->mu = mu;

  (void)vtrust_state;

  return GSL_SUCCESS;
}

static int
svd_solve(const gsl_vector * f, gsl_vector *x,
          const void * vtrust_state, void *vstate)
{
  int status = GSL_SUCCESS;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  svd_state_t *state = (svd_state_t *) vstate;
  const size_t p = state->p;
  const double tol = GSL_DBL_EPSILON;
  const double s0 = gsl_vector_get(state->S, 0);
  size_t j;

  /* compute workp = - U^T f */
  gsl_blas_dgemv(CblasTrans, -1.0, state->U, f, 0.0, state->workp);

  /*
   * compute:
   *
   * workp = sum_i s_i / (s_i^2 + mu) (-u_i^T f)
   */

  if (state->mu == 0.0)
    {
      /*
       * compute Gauss-Newton direction by solving
       * J x = -f
       */

      for (j = 0; j < p; ++j)
        {
          double sj = gsl_vector_get(state->S, j);
          double *ptr = gsl_vector_ptr(state->workp, j);
          double alpha;

          if (sj <= tol * s0)
            alpha = 0.0;
          else
            alpha = 1.0 / sj;

          *ptr *= alpha;
        }
    }
  else
    {
      /*
       * solve:
       *
       * [  J D^{-1}  ] (D x) = -[ f ]
       * [ sqrt(mu) I ]          [ 0 ]
       *
       * using SVD factorization of J D^{-1}
       */

      for (j = 0; j < p; ++j)
        {
          double sj = gsl_vector_get(state->S, j);
          double *ptr = gsl_vector_ptr(state->workp, j);

          *ptr *= sj / (sj*sj + state->mu);
        }
    }

  /* compute: x = V * workp */
  gsl_blas_dgemv(CblasNoTrans, 1.0, state->V, state->workp, 0.0, x);

  /* compute D^{-1} x */
  gsl_vector_div(x, trust_state->diag);

  return status;
}

static int
svd_rcond(double * rcond, void * vstate)
{
  int status = GSL_SUCCESS;
  svd_state_t *state = (svd_state_t *) vstate;
  double smax = gsl_vector_get(state->S, 0);
  double smin = gsl_vector_get(state->S, state->p - 1);

  *rcond = smin / smax;

  return status;
}

static const gsl_multifit_nlinear_solver svd_type =
{
  "svd",
  svd_alloc,
  svd_init,
  svd_presolve,
  svd_solve,
  svd_rcond,
  svd_free
};

const gsl_multifit_nlinear_solver *gsl_multifit_nlinear_solver_svd = &svd_type;
