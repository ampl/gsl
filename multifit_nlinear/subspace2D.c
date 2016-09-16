/* multifit_nlinear/subspace2D.c
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
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>

#define SCALE_SUB2D     1

/*
 * This module implements a 2D subspace trust region subproblem method,
 * as outlined in
 *
 * [1] G. A. Shultz, R. B. Schnabel, and R. H. Byrd
 *     A Family of Trust-Region-Based Algorithms for Unconstrained
 *     Minimization with Strong Global Convergence Properties,
 *     SIAM Journal on Numerical Analysis 1985 22:1, 47-67 
 *
 * [2] R. H. Byrd, R. B. Schnabel, G. A. Shultz,
 *     Approximate solution of the trust region problem by
 *     minimization over two-dimensional subspaces,
 *     Mathematical Programming, January 1988, Volume 40,
 *     Issue 1, pp 247-263
 *
 * The idea is to solve:
 *
 * min_{dx} g^T dx + 1/2 dx^T B dx
 * 
 * with constraints:
 *
 * ||dx|| <= delta
 * dx \in span{dx_sd, dx_gn}
 *
 * where B is the Hessian matrix, B = J^T J
 *
 * The steps are as follows:
 *
 * 1. preloop:
 *    a. Compute Gauss-Newton and steepest descent vectors,
 *       dx_gn, dx_sd
 *    b. Compute an orthonormal basis for span(dx_sd, dx_gn) by
 *       constructing W = [ dx_sd, dx_gn ] and performing a QR
 *       decomposition of W. The 2 columns of the Q matrix
 *       will then span the column space of W. W should have rank 2
 *       unless dx_sd and dx_gn are parallel, in which case it will
 *       have rank 1.
 *    c. Precompute various quantities needed for the step calculation
 *
 * 2. step:
 *    a. If the Gauss-Newton step is inside the trust region, use it
 *    b. if W has rank 1, we cannot form a 2D subspace, so in this case
 *       follow the steepest descent direction to the trust region boundary
 *       and use that as the step.
 *    c. In the full rank 2 case, if the GN point is outside the trust region,
 *       then the minimizer of the objective function lies on the trust
 *       region boundary. Therefore the minimization problem becomes:
 *
 *       min_{dx} g^T dx + 1/2 dx^T B dx, with ||dx|| = delta, dx = Q * x
 *
 *       where x is a 2-vector to be determined and the columns of Q are
 *       the orthonormal basis vectors of the subspace. Note the equality
 *       constraint now instead of <=. In terms of the new variable x,
 *       the minimization problem becomes:
 *
 *       min_x subg^T x + 1/2 x^T subB x, with ||Q*x|| = ||x|| = delta
 *
 *       where:
 *         subg = Q^T g   (2-by-1)
 *         subB = Q^T B Q (2-by-2)
 *
 *       This equality constrained 2D minimization problem can be solved
 *       with a Lagrangian multiplier, which results in a 4th degree polynomial
 *       equation to be solved. The equation is:
 *
 *         lambda^4  1
 *       + lambda^3  2 tr(B)
 *       + lambda^2  (tr(B)^2 + 2 det(B) - g^T g / delta^2)
 *       + lambda^1  (2 det(B) tr(B) - 2 g^T adj(B)^T g / delta^2)
 *       + lambda^0  (det(B)^2 - g^T adj(B)^T adj(B) g / delta^2)
 *
 *       where adj(B) is the adjugate matrix of B.
 *
 *       We then check each of the 4 solutions for lambda to determine which
 *       lambda results in the smallest objective function value. This x
 *       is then used to construct the final step: dx = Q*x
 */

typedef struct
{
  size_t n;                  /* number of observations */
  size_t p;                  /* number of parameters */
  gsl_vector *dx_gn;         /* Gauss-Newton step, size p */
  gsl_vector *dx_sd;         /* steepest descent step, size p */
  double norm_Dgn;           /* || D dx_gn || */
  double norm_Dsd;           /* || D dx_sd || */
  gsl_vector *workp;         /* workspace, length p */
  gsl_vector *workn;         /* workspace, length n */
  gsl_matrix *W;             /* orthonormal basis for 2D subspace, p-by-2 */
  gsl_matrix *JQ;            /* J * Q, n-by-p */
  gsl_vector *tau;           /* Householder scalars */
  gsl_vector *subg;          /* subspace gradient = W^T g, 2-by-1 */
  gsl_matrix *subB;          /* subspace Hessian = W^T B W, 2-by-2 */
  gsl_permutation *perm;     /* permutation matrix */

  double trB;                /* Tr(subB) */
  double detB;               /* det(subB) */
  double normg;              /* || subg || */
  double term0;              /* g^T adj(B)^T adj(B) g */
  double term1;              /* g^T adj(B)^T g */

  size_t rank;               /* rank of [ dx_sd, dx_gn ] matrix */

  gsl_poly_complex_workspace *poly_p;

  /* tunable parameters */
  gsl_multifit_nlinear_parameters params;
} subspace2D_state_t;

#include "common.c"

static void * subspace2D_alloc (const void * params, const size_t n, const size_t p);
static void subspace2D_free(void *vstate);
static int subspace2D_init(const void *vtrust_state, void *vstate);
static int subspace2D_preloop(const void * vtrust_state, void * vstate);
static int subspace2D_step(const void * vtrust_state, const double delta,
                           gsl_vector * dx, void * vstate);
static int subspace2D_preduction(const void * vtrust_state, const gsl_vector * dx,
                                 double * pred, void * vstate);
static int subspace2D_solution(const double lambda, gsl_vector * x,
                               subspace2D_state_t * state);
static double subspace2D_objective(const gsl_vector * x, subspace2D_state_t * state);

static void *
subspace2D_alloc (const void * params, const size_t n, const size_t p)
{
  const gsl_multifit_nlinear_parameters *par = (const gsl_multifit_nlinear_parameters *) params;
  subspace2D_state_t *state;
  
  state = calloc(1, sizeof(subspace2D_state_t));
  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate subspace2D state", GSL_ENOMEM);
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

  state->W = gsl_matrix_alloc(p, 2);
  if (state->W == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for W", GSL_ENOMEM);
    }

  state->JQ = gsl_matrix_alloc(n, p);
  if (state->JQ == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for JQ", GSL_ENOMEM);
    }

  state->tau = gsl_vector_alloc(2);
  if (state->tau == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for tau", GSL_ENOMEM);
    }

  state->subg = gsl_vector_alloc(2);
  if (state->subg == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for subg", GSL_ENOMEM);
    }

  state->subB = gsl_matrix_alloc(2, 2);
  if (state->subB == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for subB", GSL_ENOMEM);
    }

  state->perm = gsl_permutation_alloc(2);
  if (state->perm == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for perm", GSL_ENOMEM);
    }

  state->poly_p = gsl_poly_complex_workspace_alloc(5);
  if (state->poly_p == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for poly workspace", GSL_ENOMEM);
    }

  state->n = n;
  state->p = p;
  state->rank = 0;
  state->params = *par;

  return state;
}

static void
subspace2D_free(void *vstate)
{
  subspace2D_state_t *state = (subspace2D_state_t *) vstate;

  if (state->dx_gn)
    gsl_vector_free(state->dx_gn);

  if (state->dx_sd)
    gsl_vector_free(state->dx_sd);

  if (state->workp)
    gsl_vector_free(state->workp);

  if (state->workn)
    gsl_vector_free(state->workn);

  if (state->W)
    gsl_matrix_free(state->W);

  if (state->JQ)
    gsl_matrix_free(state->JQ);

  if (state->tau)
    gsl_vector_free(state->tau);

  if (state->subg)
    gsl_vector_free(state->subg);

  if (state->subB)
    gsl_matrix_free(state->subB);

  if (state->perm)
    gsl_permutation_free(state->perm);

  if (state->poly_p)
    gsl_poly_complex_workspace_free(state->poly_p);

  free(state);
}

/*
subspace2D_init()
  Initialize subspace2D solver

Inputs: vtrust_state - trust state
        vstate       - workspace

Return: success/error
*/

static int
subspace2D_init(const void *vtrust_state, void *vstate)
{
  (void)vtrust_state;
  (void)vstate;

  return GSL_SUCCESS;
}

/*
subspace2D_preloop()
  Initialize subspace2D method prior to iteration loop.
This involves computing the Gauss-Newton step and
steepest descent step

Notes: on output,
1) state->dx_gn contains Gauss-Newton step
2) state->dx_sd contains steepest descent step
3) state->rank contains the rank([dx_sd, dx_gn])
4) if full rank subspace (rank = 2), then:
   state->trB = Tr(subB)
   state->detB = det(subB)
   state->normg = || subg ||
*/

static int
subspace2D_preloop(const void * vtrust_state, void * vstate)
{
  int status;
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  subspace2D_state_t *state = (subspace2D_state_t *) vstate;
  const gsl_multifit_nlinear_parameters *params = trust_state->params;
  double u;
  double norm_g;  /* ||g|| */
  double norm_Jg; /* || J g || */
  double alpha;   /* ||g||^2 / ||Jg||^2 */
  double norm_Dinvg;   /* || D^{-1} g || */
  double norm_JDinv2g; /* || J D^{-2} g || */
  gsl_vector_view v;
  double work_data[2];
  gsl_vector_view work = gsl_vector_view_array(work_data, 2);
  int signum;

#if 0 /* XXX */
  print_octave(trust_state->J, "J");
  printv_octave(trust_state->f, "f");
  printv_octave(trust_state->diag, "d");
#endif

  /* initialize linear least squares solver */
  status = (params->solver->init)(trust_state, trust_state->solver_state);
  if (status)
    return status;

  /* prepare the linear solver to compute Gauss-Newton step */
  status = (params->solver->presolve)(0.0, trust_state, trust_state->solver_state);
  if (status)
    return status;

  /* solve: J dx_gn = -f for Gauss-Newton step */
  status = (params->solver->solve)(trust_state->f,
                                   state->dx_gn,
                                   trust_state,
                                   trust_state->solver_state);
  if (status)
    return status;

  /* now calculate the steepest descent step */

#if SCALE_SUB2D

  /* compute: workn = J*g */
  gsl_blas_dgemv(CblasNoTrans, 1.0, trust_state->J, trust_state->g, 0.0, state->workn);

  /* compute |g| and |Jg| */
  norm_g = gsl_blas_dnrm2(trust_state->g);
  norm_Jg = gsl_blas_dnrm2(state->workn);

  /* alpha = |g|^2 / |Jg|^2 */
  u = norm_g / norm_Jg;
  alpha = u * u;

  /* dx_sd = -alpha * g */
  gsl_vector_memcpy(state->dx_sd, trust_state->g);
  gsl_vector_scale(state->dx_sd, -alpha);

  /* store norms */
  state->norm_Dgn = gsl_blas_dnrm2(state->dx_gn);
  state->norm_Dsd = gsl_blas_dnrm2(state->dx_sd);

#else

  /* compute workp = D^{-1} g and its norm */
  gsl_vector_memcpy(state->workp, trust_state->g);
  gsl_vector_div(state->workp, trust_state->diag);
  norm_Dinvg = gsl_blas_dnrm2(state->workp);

  /* compute workp = D^{-2} g */
  gsl_vector_div(state->workp, trust_state->diag);

  /* compute: workn = J D^{-2} g */
  gsl_blas_dgemv(CblasNoTrans, 1.0, trust_state->J, state->workp, 0.0, state->workn);
  norm_JDinv2g = gsl_blas_dnrm2(state->workn);

  u = norm_Dinvg / norm_JDinv2g;
  alpha = u * u;

  /* dx_sd = -alpha D^{-2} g */
  gsl_vector_memcpy(state->dx_sd, state->workp);
  gsl_vector_scale(state->dx_sd, -alpha);

  /* store norms */
  state->norm_Dgn = scaled_enorm(trust_state->diag, state->dx_gn);
  state->norm_Dsd = scaled_enorm(trust_state->diag, state->dx_sd);

#endif

  /*
   * now compute orthonormal basis for span(dx_sd, dx_gn) using
   * QR decomposition; set W = [ dx_sd, dx_gn ] and then
   * the Q matrix will form a basis for Col(W)
   */

  v = gsl_matrix_column(state->W, 0);
  gsl_vector_memcpy(&v.vector, state->dx_sd);
  u = gsl_blas_dnrm2(state->dx_sd);
  if (u != 0)
    gsl_vector_scale(&v.vector, 1.0 / u);

  v = gsl_matrix_column(state->W, 1);
  gsl_vector_memcpy(&v.vector, state->dx_gn);
  u = gsl_blas_dnrm2(state->dx_gn);
  if (u != 0)
    gsl_vector_scale(&v.vector, 1.0 / u);

#if 0
  print_octave(state->W, "W");
#endif

  /* use a rank revealing QR decomposition in case dx_sd and dx_gn
   * are parallel */
  gsl_linalg_QRPT_decomp(state->W, state->tau, state->perm, &signum, &work.vector);

  /* check for parallel dx_sd, dx_gn, in which case rank will be 1 */
#if 1 /*XXX*/
  state->rank = gsl_linalg_QRPT_rank(state->W, -1.0);
#else
  state->rank = gsl_linalg_QRPT_rank(state->W, GSL_DBL_EPSILON);
#endif

  if (state->rank == 2)
    {
      /*
       * full rank subspace, compute:
       * subg = Q^T D^{-1} g
       * subB = Q^T D^{-1} B D^{-1} Q where B = J^T J
       */
      const size_t p = state->p;
      size_t i;
      gsl_matrix_view JQ = gsl_matrix_submatrix(state->JQ, 0, 0, state->n, GSL_MIN(2, p));
      double B00, B10, B11, g0, g1;

#if SCALE_SUB2D

      /* compute subg */

      gsl_vector_memcpy(state->workp, trust_state->g);
      gsl_linalg_QR_QTvec(state->W, state->tau, state->workp);

      for (i = 0; i < 2; ++i)
        {
          double gi = gsl_vector_get(state->workp, i);
          gsl_vector_set(state->subg, i, gi);
        }

      /* JQ = J * Q */
      gsl_matrix_memcpy(state->JQ, trust_state->J);
      gsl_linalg_QR_matQ(state->W, state->tau, state->JQ);

      /* subB = (JQ)' (JQ) = Q' J' J Q = Q' B Q */
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &JQ.matrix, 0.0, state->subB);

      B00 = gsl_matrix_get(state->subB, 0, 0);
      B10 = gsl_matrix_get(state->subB, 1, 0);
      B11 = gsl_matrix_get(state->subB, 1, 1);

      g0 = gsl_vector_get(state->subg, 0);
      g1 = gsl_vector_get(state->subg, 1);

      state->trB = B00 + B11;
      state->detB = B00*B11 - B10*B10;
      state->normg = gsl_blas_dnrm2(state->subg);

      /* g^T adj(B)^T adj(B) g */
      state->term0 = (B10*B10 + B11*B11)*g0*g0 -
                     2*B10*(B00 + B11)*g0*g1 +
                     (B00*B00 + B10*B10)*g1*g1;

      /* g^T adj(B)^T g */
      state->term1 = B11 * g0 * g0 + g1 * (B00*g1 - 2*B10*g0);

#else

      /* compute subg */
      gsl_vector_memcpy(state->workp, trust_state->g);
      gsl_vector_div(state->workp, trust_state->diag);
      gsl_linalg_QR_QTvec(state->W, state->tau, state->workp);

      g0 = gsl_vector_get(state->workp, 0);
      g1 = gsl_vector_get(state->workp, 1);

      gsl_vector_set(state->subg, 0, g0);
      gsl_vector_set(state->subg, 1, g1);

      /* compute subB */

      /* compute J D^{-1} */
      gsl_matrix_memcpy(state->JQ, trust_state->J);

      for (i = 0; i < p; ++i)
        {
          gsl_vector_view c = gsl_matrix_column(state->JQ, i);
          double di = gsl_vector_get(trust_state->diag, i);
          gsl_vector_scale(&c.vector, 1.0 / di);
        }

      /* compute J D^{-1} Q */
      gsl_linalg_QR_matQ(state->W, state->tau, state->JQ);

      /* compute subB = Q^T D^{-1} J^T J D^{-1} Q */
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &JQ.matrix, 0.0, state->subB);

      B00 = gsl_matrix_get(state->subB, 0, 0);
      B10 = gsl_matrix_get(state->subB, 1, 0);
      B11 = gsl_matrix_get(state->subB, 1, 1);

      state->trB = B00 + B11;
      state->detB = B00*B11 - B10*B10;
      state->normg = gsl_blas_dnrm2(state->subg);

      /* g^T adj(B)^T adj(B) g */
      state->term0 = (B10*B10 + B11*B11)*g0*g0 -
                     2*B10*(B00 + B11)*g0*g1 +
                     (B00*B00 + B10*B10)*g1*g1;

      /* g^T adj(B)^T g */
      state->term1 = B11 * g0 * g0 + g1 * (B00*g1 - 2*B10*g0);

#endif

    }

#if 0 /*XXX*/
  fprintf(stdout, "|dx_gn| = %.12e |dx_sd| = %.12e\n", state->norm_Dgn, state->norm_Dsd);
#endif

#if 0
  printv_octave(state->dx_gn, "dx_gn");
  printv_octave(state->dx_sd, "dx_sd");
  print_octave(state->subB, "B");
  printv_octave(state->subg, "g");
#endif

  return GSL_SUCCESS;
}


/*
2Dsubspace_step()
  Calculate a new step with 2D subspace method. Based on [1]. We
seek a vector dx in span{dx_gn, dx_sd} which minimizes the model
function subject to ||dx|| <= delta
*/

static int
subspace2D_step(const void * vtrust_state, const double delta,
                gsl_vector * dx, void * vstate)
{
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  subspace2D_state_t *state = (subspace2D_state_t *) vstate;

  if (state->norm_Dgn <= delta)
    {
      /* Gauss-Newton step is inside trust region, use it as final step
       * since it is the global minimizer of the quadratic model function */
      gsl_vector_memcpy(dx, state->dx_gn);
    }
  else if (state->rank < 2)
    {
      /* rank of [dx_sd, dx_gn] is 1, meaning dx_sd and dx_gn
       * are parallel so we can't form a 2D subspace. Follow the steepest
       * descent direction to the trust region boundary as our step */
      gsl_vector_memcpy(dx, state->dx_sd);
      gsl_vector_scale(dx, delta / state->norm_Dsd);
    }
  else
    {
      int status;
      const double delta_sq = delta * delta;
      double TrB_D = state->trB * delta;
      double detB_D = state->detB * delta;
      double normg_sq = state->normg * state->normg;
      double u = state->normg / delta;
      double a[5];
      double z[8];

#if 0
      a[0] = state->detB * state->detB - state->term0 / delta_sq;
      a[1] = 2 * state->detB * state->trB - 2 * state->term1 / delta_sq;
      a[2] = state->trB * state->trB + 2 * state->detB - u * u;
      a[3] = 2 * state->trB;
      a[4] = 1.0;
#else
      a[0] = detB_D * detB_D - state->term0;
      a[1] = 2 * state->detB * state->trB * delta_sq - 2 * state->term1;
      a[2] = TrB_D * TrB_D + 2 * state->detB * delta_sq - normg_sq;
      a[3] = 2 * state->trB * delta_sq;
      a[4] = delta_sq;
#endif

      status = gsl_poly_complex_solve(a, 5, state->poly_p, z);
      if (status == GSL_SUCCESS)
        {
          size_t i;
          double min = 0.0;
          int mini = -1;
          double x_data[2];
          gsl_vector_view x = gsl_vector_view_array(x_data, 2);

          /*
           * loop through all four values of the Lagrange multiplier
           * lambda. For each lambda, evaluate the objective function
           * with Re(lambda) to determine which lambda minimizes the
           * function
           */
          for (i = 0; i < 4; ++i)
            {
              double cost, normx;

              /*fprintf(stderr, "root: %.12e + %.12e i\n",
                      z[2*i], z[2*i+1]);*/

              status = subspace2D_solution(z[2*i], &x.vector, state);
              if (status != GSL_SUCCESS)
                continue; /* singular matrix system */

              /* ensure ||x|| = delta */
              normx = gsl_blas_dnrm2(&x.vector);
              if (normx == 0.0)
                continue;

              gsl_vector_scale(&x.vector, delta / normx);

              /* evaluate objective function to determine minimizer */
              cost = subspace2D_objective(&x.vector, state);
              if (mini < 0 || cost < min)
                {
                  mini = (int) i;
                  min = cost;
                }

              /*fprintf(stderr, "|x| = %.12e (%.12e) cost = %.12e\n", gsl_blas_dnrm2(&x.vector), delta, cost);*/
            }

          if (mini < 0)
            {
              /*XXXfprintf(stderr, "ERROR: did not find minimizer\n");*/
              return GSL_FAILURE;
            }
          else
            {
              /* compute x which minimizes objective function */
              subspace2D_solution(z[2*mini], &x.vector, state);

              /* dx = W * x */
              gsl_vector_set_zero(dx);
              gsl_vector_set(dx, 0, gsl_vector_get(&x.vector, 0));
              gsl_vector_set(dx, 1, gsl_vector_get(&x.vector, 1));
              gsl_linalg_QR_Qvec(state->W, state->tau, dx);

#if !SCALE_SUB2D
              gsl_vector_div(dx, trust_state->diag);
#endif
            }
        }
      else
        {
          GSL_ERROR ("gsl_poly_complex_solve failed", status);
        }
    }

#if 0
#if SCALE_SUB2D
  fprintf(stdout, "|dx| = %.12e\n", gsl_blas_dnrm2(dx));
#else
  fprintf(stdout, "|dx| = %.12e\n", scaled_enorm(trust_state->diag, dx));
#endif
#endif

  return GSL_SUCCESS;
}

static int
subspace2D_preduction(const void * vtrust_state, const gsl_vector * dx,
                  double * pred, void * vstate)
{
  const gsl_multifit_nlinear_trust_state *trust_state =
    (const gsl_multifit_nlinear_trust_state *) vtrust_state;
  subspace2D_state_t *state = (subspace2D_state_t *) vstate;

  *pred = quadratic_preduction(trust_state->f, trust_state->J, dx, state->workn);

  return GSL_SUCCESS;
}

/* solve 2D subspace problem: (B + lambda*I) x = -g */
static int
subspace2D_solution(const double lambda, gsl_vector * x,
                    subspace2D_state_t * state)
{
  int status = GSL_SUCCESS;
  double C_data[4];
  gsl_matrix_view C = gsl_matrix_view_array(C_data, 2, 2);
  double B00 = gsl_matrix_get(state->subB, 0, 0);
  double B10 = gsl_matrix_get(state->subB, 1, 0);
  double B11 = gsl_matrix_get(state->subB, 1, 1);
  int signum;

  /* construct C = B + lambda*I */
  gsl_matrix_set(&C.matrix, 0, 0, B00 + lambda);
  gsl_matrix_set(&C.matrix, 1, 0, B10);
  gsl_matrix_set(&C.matrix, 0, 1, B10);
  gsl_matrix_set(&C.matrix, 1, 1, B11 + lambda);

#if 0
  gsl_linalg_mcholesky_decomp(&C.matrix, state->perm, NULL);
  gsl_linalg_mcholesky_solve(&C.matrix, state->perm, state->subg, x);
#elif 1
  {
    double D_data[4];
    double V_data[4];
    double S_data[2];
    double work_data[2];
    gsl_matrix_view V = gsl_matrix_view_array(V_data, 2, 2);
    gsl_matrix_view D = gsl_matrix_view_array(D_data, 2, 2);
    gsl_vector_view S = gsl_vector_view_array(S_data, 2);
    gsl_vector_view work = gsl_vector_view_array(work_data, 2);
    double tol;
    double S0, S1;
    double r;

    gsl_matrix_memcpy(&D.matrix, &C.matrix);

#if 0
    gsl_linalg_SV_decomp(&C.matrix, &V.matrix, &S.vector, &work.vector);
#else
    gsl_linalg_SV_decomp_jacobi(&C.matrix, &V.matrix, &S.vector);
#endif

    S0 = S_data[0];
    S1 = S_data[1];
    tol = 2 * S0 * GSL_DBL_EPSILON;

    if (S1 < tol)
      {
        /* matrix is near-singular, don't try to compute solution */
        return GSL_ESING;
      }

    gsl_linalg_SV_solve(&C.matrix, &V.matrix, &S.vector, state->subg, x);

    gsl_vector_memcpy(&work.vector, state->subg);
    gsl_blas_dsymv(CblasUpper, -1.0, &D.matrix, x, 1.0, &work.vector);

    r = gsl_blas_dnrm2(&work.vector);

#if 0 /* XXX */
    if (S1 < tol && r > 1.0e3)
      {
        fprintf(stderr, "residual = %.12e\n", r);
        print_octave(&D.matrix, "C");
        printv_octave(state->subg, "g");
        printv_octave(x, "y");

        gsl_linalg_mcholesky_decomp(&C.matrix, state->perm, NULL);
        gsl_linalg_mcholesky_solve(&C.matrix, state->perm, state->subg, x);
        printv_octave(x, "y_chol");

        exit(1);
      }
#endif
  }
#endif

  gsl_vector_scale(x, -1.0);

  return status;
}

/* evaluate 2D objective function: f(x) = g^T x + 1/2 x^T B x */
static double
subspace2D_objective(const gsl_vector * x, subspace2D_state_t * state)
{
  double f;
  double y_data[2];
  gsl_vector_view y = gsl_vector_view_array(y_data, 2);

  /* compute: y = g + 1/2 B x */
  gsl_vector_memcpy(&y.vector, state->subg);
  gsl_blas_dsymv(CblasLower, 0.5, state->subB, x, 1.0, &y.vector);

  /* compute: f = x^T ( g + 1/2 B x ) */
  gsl_blas_ddot(x, &y.vector, &f);

  return f;
}

static const gsl_multifit_nlinear_trs subspace2D_type =
{
  "2D-subspace",
  subspace2D_alloc,
  subspace2D_init,
  subspace2D_preloop,
  subspace2D_step,
  subspace2D_preduction,
  subspace2D_free
};

const gsl_multifit_nlinear_trs *gsl_multifit_nlinear_trs_subspace2D = &subspace2D_type;
