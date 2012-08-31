/* ode-initval2/modnewton1.c
 * 
 * Copyright (C) 2008, 2009, 2010 Tuomo Keskitalo
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

/* A modified Newton iteration method for solving non-linear 
   equations in implicit Runge-Kutta methods.
*/

/* References: 

   Ascher, U.M., Petzold, L.R., Computer methods for ordinary
   differential and differential-algebraic equations, SIAM, 
   Philadelphia, 1998. ISBN 0898714125, 9780898714128
   
   Hairer, E., Wanner, G., Solving Ordinary Differential 
   Equations II: Stiff and Differential-Algebraic Problems,
   Springer, 1996. ISBN 3540604529, 9783540604525
*/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "odeiv_util.h"

typedef struct
{
  /* iteration matrix I - h A (*) J */
  gsl_matrix *IhAJ;

  /* permutation for LU-decomposition */
  gsl_permutation *p;

  /* difference vector for kth Newton iteration */
  gsl_vector *dYk;

  /* scaled dYk with desired error level */
  gsl_vector *dScal;

  /* current Runge-Kutta points (absolute values) */
  double *Yk;

  /* function values at points Yk */
  double *fYk;

  /* right hand side of Newton iteration formula */
  gsl_vector *rhs;

  /* stopping criterion value from previous step */
  double eeta_prev;
}
modnewton1_state_t;

static void *
modnewton1_alloc (size_t dim, size_t stage)
{
  modnewton1_state_t *state =
    (modnewton1_state_t *) malloc (sizeof (modnewton1_state_t));

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for modnewton1_state_t",
                      GSL_ENOMEM);
    }

  state->IhAJ = gsl_matrix_alloc (dim * stage, dim * stage);

  if (state->IhAJ == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for IhAJ", GSL_ENOMEM);
    }

  state->p = gsl_permutation_alloc (dim * stage);

  if (state->p == 0)
    {
      gsl_matrix_free (state->IhAJ);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for p", GSL_ENOMEM);
    }

  state->dYk = gsl_vector_alloc (dim * stage);

  if (state->dYk == 0)
    {
      gsl_permutation_free (state->p);
      gsl_matrix_free (state->IhAJ);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for dYk", GSL_ENOMEM);
    }

  state->dScal = gsl_vector_alloc (dim * stage);

  if (state->dScal == 0)
    {
      gsl_vector_free (state->dYk);
      gsl_permutation_free (state->p);
      gsl_matrix_free (state->IhAJ);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for dScal", GSL_ENOMEM);
    }

  state->Yk = (double *) malloc (dim * stage * sizeof (double));

  if (state->Yk == 0)
    {
      gsl_vector_free (state->dScal);
      gsl_vector_free (state->dYk);
      gsl_permutation_free (state->p);
      gsl_matrix_free (state->IhAJ);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for Yk", GSL_ENOMEM);
    }

  state->fYk = (double *) malloc (dim * stage * sizeof (double));

  if (state->fYk == 0)
    {
      free (state->Yk);
      gsl_vector_free (state->dScal);
      gsl_vector_free (state->dYk);
      gsl_permutation_free (state->p);
      gsl_matrix_free (state->IhAJ);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for Yk", GSL_ENOMEM);
    }

  state->rhs = gsl_vector_alloc (dim * stage);

  if (state->rhs == 0)
    {
      free (state->fYk);
      free (state->Yk);
      gsl_vector_free (state->dScal);
      gsl_vector_free (state->dYk);
      gsl_permutation_free (state->p);
      gsl_matrix_free (state->IhAJ);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for rhs", GSL_ENOMEM);
    }

  state->eeta_prev = GSL_DBL_MAX;

  return state;
}

static int
modnewton1_init (void *vstate, const gsl_matrix * A,
                 const double h, const gsl_matrix * dfdy,
                 const gsl_odeiv2_system * sys)
{
  /* Initializes the method by forming the iteration matrix IhAJ
     and generating its LU-decomposition
   */

  modnewton1_state_t *state = (modnewton1_state_t *) vstate;

  gsl_matrix *const IhAJ = state->IhAJ;
  gsl_permutation *const p = state->p;

  const size_t dim = sys->dimension;
  const size_t stage = A->size1;

  state->eeta_prev = GSL_DBL_MAX;

  /* Generate IhAJ */

  {
    size_t i, j, k, m;

    for (i = 0; i < dim; i++)
      for (j = 0; j < dim; j++)
        for (k = 0; k < stage; k++)
          for (m = 0; m < stage; m++)
            {
              size_t x = dim * k + i;
              size_t y = dim * m + j;

              if (x != y)
                gsl_matrix_set (IhAJ, x, y,
                                -h * gsl_matrix_get (A, k, m) *
                                gsl_matrix_get (dfdy, i, j));
              else
                gsl_matrix_set (IhAJ, x, y,
                                1.0 - h * gsl_matrix_get (A, k, m) *
                                gsl_matrix_get (dfdy, i, j));
            }
  }

  /* decompose */

  {
    int signum;
    int s = gsl_linalg_LU_decomp (IhAJ, p, &signum);

    if (s != GSL_SUCCESS)
      return s;
  }

  return GSL_SUCCESS;
}

static int
modnewton1_solve (void *vstate, const gsl_matrix * A,
                  const double c[], const double t, const double h,
                  const double y0[], const gsl_odeiv2_system * sys,
                  double YZ[], const double errlev[])
{
  /* Solves the non-linear equation system resulting from implicit
     Runge-Kutta methods by a modified Newton iteration. The
     modified Newton iteration with this problem is stated in the
     form

     IhAJ * dYk = rhs

     in which IhAJ is the iteration matrix: IhAJ = (I - hA (*) J) in
     which (*) is the Kronecker matrix product (size of IhAJ =
     dim*stage, dim*stage), dYk is the Runge-Kutta point (Y)
     difference vector for kth Newton iteration: dYk = Y(k+1) - Y(k),
     and rhs = Y(k) - y0 - h * sum j=1..stage (a_j * f(Y(k)))

     This function solves dYk by LU-decomposition of IhAJ with partial
     pivoting.
   */

  modnewton1_state_t *state = (modnewton1_state_t *) vstate;

  gsl_matrix *const IhAJ = state->IhAJ;
  gsl_permutation *const p = state->p;
  gsl_vector *const dYk = state->dYk;
  double *const Yk = state->Yk;
  double *const fYk = state->fYk;
  gsl_vector *const rhs = state->rhs;
  double *const eeta_prev = &(state->eeta_prev);
  gsl_vector *const dScal = state->dScal;

  const size_t dim = sys->dimension;
  const size_t stage = A->size1;

  int status = GSL_CONTINUE;    /* Convergence status for Newton iteration */

  /* Set starting values for iteration. Use starting values of Y(k) =
     y0. FIXME: Implement better initial guess described in Hairer and
     Wanner.
   */

  {
    size_t j, m;

    gsl_vector_set_zero (dYk);

    for (j = 0; j < stage; j++)
      for (m = 0; m < dim; m++)
        Yk[j * dim + m] = y0[m];
  }

  /* Loop modified Newton iterations */

  {
    size_t iter = 0;
    size_t j, m, n;

    /* Maximum number of Newton iterations. */
    const size_t max_iter = 7;

    double dScal_norm = 0.0;    /* Newton iteration step length */
    double dScal_norm_prev = 0.0;       /* Previous dScal_norm */

    while (status == GSL_CONTINUE && iter < max_iter)
      {
        iter++;

        /* Update Y(k) and evaluate function */

        for (j = 0; j < stage; j++)
          {
            for (m = 0; m < dim; m++)
              {
                Yk[j * dim + m] += gsl_vector_get (dYk, j * dim + m);
              }

            {
              int s = GSL_ODEIV_FN_EVAL (sys, t + c[j] * h, &Yk[j * dim],
                                         &fYk[j * dim]);
              if (s != GSL_SUCCESS)
                {
                  return s;
                }
            }
          }

        /* Calculate rhs  */

        for (j = 0; j < stage; j++)
          for (m = 0; m < dim; m++)
            {
              double sum = 0;

              for (n = 0; n < stage; n++)
                sum += gsl_matrix_get (A, j, n) * fYk[n * dim + m];

              gsl_vector_set (rhs, j * dim + m,
                              -1.0 * Yk[j * dim + m] + y0[m] + h * sum);
            }

        /* Solve dYk */

        {
          int s = gsl_linalg_LU_solve (IhAJ, p, rhs, dYk);

          if (s != GSL_SUCCESS)
            {
              return s;
            }
        }

        /* Evaluate convergence according to method by Hairer and
           Wanner, section IV.8. The iteration step is normalized by
           desired error level before calculation of the norm, to take
           the tolerance on each component of y into account.
         */

        {
          /* relative change in two Newton iteration steps */
          double theta_k;

          /* transformation of theta_k */
          double eeta_k = 0.0;

          for (j = 0; j < stage; j++)
            for (m = 0; m < dim; m++)
              {
                gsl_vector_set (dScal, j * dim + m,
                                gsl_vector_get (dYk, j * dim + m)
                                / errlev[m]);
              }

          dScal_norm_prev = dScal_norm;
          dScal_norm = gsl_blas_dnrm2 (dScal);

          theta_k = dScal_norm / dScal_norm_prev;

          if (iter > 1)
            {
              /* check for increase in step size, which indicates divergence */

              if (theta_k >= 1.0)
                {
                  return GSL_FAILURE;
                }

              eeta_k = theta_k / (1.0 - theta_k);
            }

          else
            {
              eeta_k = pow (GSL_MAX_DBL (*eeta_prev, GSL_DBL_EPSILON), 0.8);
            }

          *eeta_prev = eeta_k;

          if (eeta_k * dScal_norm < 1.0)
            {
              status = GSL_SUCCESS;
            }
        }
      }
  }

  /* No convergence reached within allowed iterations */

  if (status != GSL_SUCCESS)
    {
      return GSL_FAILURE;
    }

  /* Give solution in YZ */

  else
    {
      size_t j, m;

      for (j = 0; j < stage; j++)
        for (m = 0; m < dim; m++)
          YZ[j * dim + m] =
            Yk[j * dim + m] + gsl_vector_get (dYk, j * dim + m);

      return status;
    }
}

static void
modnewton1_free (void *vstate)
{
  modnewton1_state_t *state = (modnewton1_state_t *) vstate;

  gsl_vector_free (state->rhs);
  free (state->fYk);
  free (state->Yk);
  gsl_vector_free (state->dScal);
  gsl_vector_free (state->dYk);
  gsl_permutation_free (state->p);
  gsl_matrix_free (state->IhAJ);
  free (state);
}
