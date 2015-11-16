/* gmres.c
 * 
 * Copyright (C) 2014 Patrick Alken
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
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_splinalg.h>

/*
 * The code in this module is based on the Householder GMRES
 * algorithm described in
 *
 * [1] H. F. Walker, Implementation of the GMRES method using
 *     Householder transformations, SIAM J. Sci. Stat. Comput.
 *     9(1), 1988.
 *
 * [2] Y. Saad, Iterative methods for sparse linear systems,
 *     2nd edition, SIAM, 2003.
 */

typedef struct
{
  size_t n;        /* size of linear system */
  size_t m;        /* dimension of Krylov subspace K_m */
  gsl_vector *r;   /* residual vector r = b - A*x */
  gsl_matrix *H;   /* Hessenberg matrix n-by-(m+1) */
  gsl_vector *tau; /* householder scalars */
  gsl_vector *y;   /* least squares rhs and solution vector */

  double *c;       /* Givens rotations */
  double *s;

  double normr;    /* residual norm ||r|| */
} gmres_state_t;

static void gmres_free(void *vstate);
static int gmres_iterate(const gsl_spmatrix *A, const gsl_vector *b,
                         const double tol, gsl_vector *x, void *vstate);

/*
gmres_alloc()
  Allocate a GMRES workspace for solving an n-by-n system A x = b

Inputs: n        - size of system
        krylov_m - size of Krylov subspace (ie: number of inner iterations)
                   if this parameter is 0, the value GSL_MIN(n,10) is
                   used

Return: pointer to workspace
*/

static void *
gmres_alloc(const size_t n, const size_t m)
{
  gmres_state_t *state;

  if (n == 0)
    {
      GSL_ERROR_NULL("matrix dimension n must be a positive integer",
                     GSL_EINVAL);
    }

  state = calloc(1, sizeof(gmres_state_t));
  if (!state)
    {
      GSL_ERROR_NULL("failed to allocate gmres state", GSL_ENOMEM);
    }

  state->n = n;

  /* compute size of Krylov subspace */
  if (m == 0)
    state->m = GSL_MIN(n, 10);
  else
    state->m = GSL_MIN(n, m);

  state->r = gsl_vector_alloc(n);
  if (!state->r)
    {
      gmres_free(state);
      GSL_ERROR_NULL("failed to allocate r vector", GSL_ENOMEM);
    }

  state->H = gsl_matrix_alloc(n, state->m + 1);
  if (!state->H)
    {
      gmres_free(state);
      GSL_ERROR_NULL("failed to allocate H matrix", GSL_ENOMEM);
    }

  state->tau = gsl_vector_alloc(state->m + 1);
  if (!state->tau)
    {
      gmres_free(state);
      GSL_ERROR_NULL("failed to allocate tau vector", GSL_ENOMEM);
    }

  state->y = gsl_vector_alloc(state->m + 1);
  if (!state->y)
    {
      gmres_free(state);
      GSL_ERROR_NULL("failed to allocate y vector", GSL_ENOMEM);
    }

  state->c = malloc(state->m * sizeof(double));
  state->s = malloc(state->m * sizeof(double));
  if (!state->c || !state->s)
    {
      gmres_free(state);
      GSL_ERROR_NULL("failed to allocate Givens vectors", GSL_ENOMEM);
    }

  state->normr = 0.0;

  return state;
} /* gmres_alloc() */

static void
gmres_free(void *vstate)
{
  gmres_state_t *state = (gmres_state_t *) vstate;

  if (state->r)
    gsl_vector_free(state->r);

  if (state->H)
    gsl_matrix_free(state->H);

  if (state->tau)
    gsl_vector_free(state->tau);

  if (state->y)
    gsl_vector_free(state->y);

  if (state->c)
    free(state->c);

  if (state->s)
    free(state->s);

  free(state);
} /* gmres_free() */

/*
gmres_iterate()
  Solve A*x = b using GMRES algorithm

Inputs: A    - sparse square matrix
        b    - right hand side vector
        tol  - stopping tolerance (see below)
        x    - (input/output) on input, initial estimate x_0;
               on output, solution vector
        work - workspace

Return:
GSL_SUCCESS if converged to solution (solution stored in x). In
this case the following will be true:

||b - A*x|| <= tol * ||b||

GSL_CONTINUE if not yet converged; in this case x contains the
most recent solution vector and calling this function more times
with the input x could result in convergence (ie: restarted GMRES)

Notes:
1) Based on algorithm 2.2 of (Walker, 1998 [1]) and algorithm 6.10 of
(Saad, 2003 [2])

2) On output, work->normr contains ||b - A*x||
*/

static int
gmres_iterate(const gsl_spmatrix *A, const gsl_vector *b,
              const double tol, gsl_vector *x,
              void *vstate)
{
  const size_t N = A->size1;
  gmres_state_t *state = (gmres_state_t *) vstate;

  if (N != A->size2)
    {
      GSL_ERROR("matrix must be square", GSL_ENOTSQR);
    }
  else if (N != b->size)
    {
      GSL_ERROR("matrix does not match right hand side", GSL_EBADLEN);
    }
  else if (N != x->size)
    {
      GSL_ERROR("matrix does not match solution vector", GSL_EBADLEN);
    }
  else if (N != state->n)
    {
      GSL_ERROR("matrix does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;
      const size_t maxit = state->m;
      const double normb = gsl_blas_dnrm2(b); /* ||b|| */
      const double reltol = tol * normb;      /* tol*||b|| */
      double normr;                           /* ||r|| */
      size_t m, k;
      double tau;                             /* householder scalar */
      gsl_matrix *H = state->H;               /* Hessenberg matrix */
      gsl_vector *r = state->r;               /* residual vector */
      gsl_vector *w = state->y;               /* least squares RHS */
      gsl_matrix_view Rm;                     /* R_m = H(1:m,2:m+1) */
      gsl_vector_view ym;                     /* y(1:m) */
      gsl_vector_view h0 = gsl_matrix_column(H, 0);

      /*
       * The Hessenberg matrix will have the following structure:
       *
       * H = [ ||r_0|| | v_1 v_2 ... v_m     ]
       *     [   u_1   | u_2 u_3 ... u_{m+1} ]
       *
       * where v_j are the orthonormal vectors spanning the Krylov
       * subpsace of length j + 1 and u_{j+1} are the householder
       * vectors of length n - j - 1.
       * In fact, u_{j+1} has length n - j since u_{j+1}[0] = 1,
       * but this 1 is not stored.
       */
      gsl_matrix_set_zero(H);

      /* Step 1a: compute r = b - A*x_0 */
      gsl_vector_memcpy(r, b);
      gsl_spblas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);

      /* Step 1b */
      gsl_vector_memcpy(&h0.vector, r);
      tau = gsl_linalg_householder_transform(&h0.vector);

      /* store tau_1 */
      gsl_vector_set(state->tau, 0, tau);

      /* initialize w (stored in state->y) */
      gsl_vector_set_zero(w);
      gsl_vector_set(w, 0, gsl_vector_get(&h0.vector, 0));

      for (m = 1; m <= maxit; ++m)
        {
          size_t j = m - 1; /* C indexing */
          double c, s;      /* Givens rotation */

          /* v_m */
          gsl_vector_view vm = gsl_matrix_column(H, m);

          /* v_m(m:end) */
          gsl_vector_view vv = gsl_vector_subvector(&vm.vector, j, N - j);

          /* householder vector u_m for projection P_m */
          gsl_vector_view um = gsl_matrix_subcolumn(H, j, j, N - j);

          /* Step 2a: form v_m = P_m e_m = e_m - tau_m w_m */
          gsl_vector_set_zero(&vm.vector);
          gsl_vector_memcpy(&vv.vector, &um.vector);
          tau = gsl_vector_get(state->tau, j); /* tau_m */
          gsl_vector_scale(&vv.vector, -tau);
          gsl_vector_set(&vv.vector, 0, 1.0 - tau);

          /* Step 2a: v_m <- P_1 P_2 ... P_{m-1} v_m */
          for (k = j; k > 0 && k--; )
            {
              gsl_vector_view uk =
                gsl_matrix_subcolumn(H, k, k, N - k);
              gsl_vector_view vk =
                gsl_vector_subvector(&vm.vector, k, N - k);
              tau = gsl_vector_get(state->tau, k);
              gsl_linalg_householder_hv(tau, &uk.vector, &vk.vector);
            }

          /* Step 2a: v_m <- A*v_m */
          gsl_spblas_dgemv(CblasNoTrans, 1.0, A, &vm.vector, 0.0, r);
          gsl_vector_memcpy(&vm.vector, r);

          /* Step 2a: v_m <- P_m ... P_1 v_m */
          for (k = 0; k <= j; ++k)
            {
              gsl_vector_view uk = gsl_matrix_subcolumn(H, k, k, N - k);
              gsl_vector_view vk = gsl_vector_subvector(&vm.vector, k, N - k);
              tau = gsl_vector_get(state->tau, k);
              gsl_linalg_householder_hv(tau, &uk.vector, &vk.vector);
            }

          /* Steps 2c,2d: find P_{m+1} and set v_m <- P_{m+1} v_m */
          if (m < N)
            {
              /* householder vector u_{m+1} for projection P_{m+1} */
              gsl_vector_view ump1 = gsl_matrix_subcolumn(H, m, m, N - m);

              tau = gsl_linalg_householder_transform(&ump1.vector);
              gsl_vector_set(state->tau, j + 1, tau);
            }

          /* Step 2e: v_m <- J_{m-1} ... J_1 v_m */
          for (k = 0; k < j; ++k)
            {
              gsl_linalg_givens_gv(&vm.vector, k, k + 1,
                                   state->c[k], state->s[k]);
            }

          if (m < N)
            {
              /* Step 2g: find givens rotation J_m for v_m(m:m+1) */
              gsl_linalg_givens(gsl_vector_get(&vm.vector, j),
                                gsl_vector_get(&vm.vector, j + 1),
                                &c, &s);

              /* store givens rotation for later use */
              state->c[j] = c;
              state->s[j] = s;

              /* Step 2h: v_m <- J_m v_m */
              gsl_linalg_givens_gv(&vm.vector, j, j + 1, c, s);

              /* Step 2h: w <- J_m w */
              gsl_linalg_givens_gv(w, j, j + 1, c, s);
            }

          /*
           * Step 2i: R_m = [ R_{m-1}, v_m ] - already taken care
           * of due to our memory storage scheme
           */

          /* Step 2j: check residual w_{m+1} for convergence */
          normr = fabs(gsl_vector_get(w, j + 1));
          if (normr <= reltol)
            {
              /*
               * method has converged, break out of loop to compute
               * update to solution vector x
               */
              break;
            }
        }

      /*
       * At this point, we have either converged to a solution or
       * completed all maxit iterations. In either case, compute
       * an update to the solution vector x and test again for
       * convergence.
       */

      /* rewind m if we exceeded maxit iterations */
      if (m > maxit)
        m--;

      /* Step 3a: solve triangular system R_m y_m = w, in place */
      Rm = gsl_matrix_submatrix(H, 0, 1, m, m);
      ym = gsl_vector_subvector(w, 0, m);
      gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit,
                     &Rm.matrix, &ym.vector);

      /*
       * Step 3b: update solution vector x; the loop below
       * uses a different but equivalent formulation from
       * Saad, algorithm 6.10, step 14; store Krylov projection
       * V_m y_m in 'r'
       */
      gsl_vector_set_zero(r);
      for (k = m; k > 0 && k--; )
        {
          double ymk = gsl_vector_get(&ym.vector, k);
          gsl_vector_view uk = gsl_matrix_subcolumn(H, k, k, N - k);
          gsl_vector_view rk = gsl_vector_subvector(r, k, N - k);

          /* r <- n_k e_k + r */
          gsl_vector_set(r, k, gsl_vector_get(r, k) + ymk);

          /* r <- P_k r */
          tau = gsl_vector_get(state->tau, k);
          gsl_linalg_householder_hv(tau, &uk.vector, &rk.vector);
        }

      /* x <- x + V_m y_m */
      gsl_vector_add(x, r);

      /* compute new residual r = b - A*x */
      gsl_vector_memcpy(r, b);
      gsl_spblas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);
      normr = gsl_blas_dnrm2(r);

      if (normr <= reltol)
        status = GSL_SUCCESS;  /* converged */
      else
        status = GSL_CONTINUE; /* not yet converged */

      /* store residual norm */
      state->normr = normr;

      return status;
    }
} /* gmres_iterate() */

static double
gmres_normr(const void *vstate)
{
  const gmres_state_t *state = (const gmres_state_t *) vstate;
  return state->normr;
} /* gmres_normr() */

static const gsl_splinalg_itersolve_type gmres_type =
{
  "gmres",
  &gmres_alloc,
  &gmres_iterate,
  &gmres_normr,
  &gmres_free
};

const gsl_splinalg_itersolve_type * gsl_splinalg_itersolve_gmres =
  &gmres_type;
