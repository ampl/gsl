/* ode-initval2/msbdf.c
 * 
 * Copyright (C) 2009, 2010 Tuomo Keskitalo
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

/* A variable-coefficient linear multistep backward differentiation
   formula (BDF) method in Nordsieck form. This stepper uses the
   explicit BDF formula as predictor and implicit BDF formula as
   corrector. A modified Newton iteration method is used to
   solve the system of non-linear equations. Method order varies
   dynamically between 1 and 5.

   References:

   Byrne, G. D., and Hindmarsh, A. C., A Polyalgorithm for the
   Numerical Solution of Ordinary Differential Equations,
   ACM Trans. Math. Software, 1 (1975), pp. 71-96.

   Brown, P. N., Byrne, G. D., and Hindmarsh, A. C., VODE: A
   Variable-coefficient ODE Solver, SIAM J. Sci. Stat. Comput. 10,
   (1989), pp. 1038-1051.

   Hindmarsh, A. C., Brown, P. N., Grant, K. E., Lee, S. L., Serban,
   R., Shumaker, D. E., and Woodward, C. S., SUNDIALS: Suite of
   Nonlinear and Differential/Algebraic Equation Solvers, ACM
   Trans. Math. Software 31 (2005), pp. 363-396.

   Note: The algorithms have been adapted for GSL ode-initval2
   framework.
*/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "odeiv_util.h"

/* Maximum order of BDF methods */
#define MSBDF_MAX_ORD 5

/* Steps until Jacobian evaluation is forced */
#define MSBDF_JAC_WAIT 50

/* Steps until iteration matrix M evaluation is forced */
#define MSBDF_M_WAIT 20

typedef struct
{
  /* Nordsieck history matrix. Includes concatenated
     Nordsieck vectors [y_n, h*y_n', (h^2/2!)*y_n'', ...,
     (h^ord/ord!)*d^(ord)(y_n)]. Nordsieck vector number i is located
     at z[i*dim] (i=0..ord).
   */
  double *z;

  double *zbackup;              /* backup of Nordsieck matrix */
  double *ytmp;                 /* work area */
  double *ytmp2;                /* work area */
  double *l;                    /* polynomial coefficients */
  double *hprev;                /* previous step sizes */
  double *hprevbackup;          /* backup of hprev */
  size_t *ordprev;              /* orders of previous calls */
  size_t *ordprevbackup;        /* backup of ordprev */
  double *errlev;               /* desired error level of y */
  gsl_vector *abscor;           /* absolute y values for correction */
  gsl_vector *abscorscaled;     /* scaled abscor for order evaluation */
  gsl_vector *relcor;           /* relative y values for correction */
  gsl_vector *svec;             /* saved abscor & work area */
  gsl_vector *tempvec;          /* work area */
  const gsl_odeiv2_driver *driver;      /* pointer to gsl_odeiv2_driver object */
  gsl_matrix *dfdy;             /* Jacobian */
  double *dfdt;                 /* storage for time derivative of f */
  gsl_matrix *M;                /* Newton iteration matrix */
  gsl_permutation *p;           /* permutation for LU decomposition of M */
  gsl_vector *rhs;              /* right hand side equations (-G) */
  long int ni;                  /* stepper call counter */
  size_t ord;                   /* current order of method */
  double tprev;                 /* t point of previous call */
  size_t ordwait;               /* counter for order change */
  size_t ordwaitbackup;         /* backup of ordwait */
  size_t failord;               /* order of convergence failure */
  double failt;                 /* t point of convergence failure */
  double ordp1coeffprev;        /* saved order coefficient */
  size_t nJ;                    /* step counter for Jacobian evaluation */
  size_t nM;                    /* step counter for update of M */
  double gammaprev;             /* gamma of previous call */
  double gammaprevbackup;       /* backup of gammaprev */
  size_t failcount;             /* counter for rejected steps */
}
msbdf_state_t;

/* Introduce msbdf_reset for use in msbdf_alloc and _apply */

static int msbdf_reset (void *, size_t);

static void *
msbdf_alloc (size_t dim)
{
  msbdf_state_t *state = (msbdf_state_t *) malloc (sizeof (msbdf_state_t));

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for msbdf_state", GSL_ENOMEM);
    }

  state->z = (double *) malloc ((MSBDF_MAX_ORD + 1) * dim * sizeof (double));

  if (state->z == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for z", GSL_ENOMEM);
    }

  state->zbackup =
    (double *) malloc ((MSBDF_MAX_ORD + 1) * dim * sizeof (double));

  if (state->zbackup == 0)
    {
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for zbackup", GSL_ENOMEM);
    }

  state->ytmp = (double *) malloc (dim * sizeof (double));

  if (state->ytmp == 0)
    {
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for ytmp", GSL_ENOMEM);
    }

  state->ytmp2 = (double *) malloc (dim * sizeof (double));

  if (state->ytmp2 == 0)
    {
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for ytmp2", GSL_ENOMEM);
    }

  state->l = (double *) malloc ((MSBDF_MAX_ORD + 1) * sizeof (double));

  if (state->l == 0)
    {
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for l", GSL_ENOMEM);
    }

  state->hprev = (double *) malloc (MSBDF_MAX_ORD * sizeof (double));

  if (state->hprev == 0)
    {
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for hprev", GSL_ENOMEM);
    }

  state->hprevbackup = (double *) malloc (MSBDF_MAX_ORD * sizeof (double));

  if (state->hprevbackup == 0)
    {
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for hprevbackup", GSL_ENOMEM);
    }

  state->ordprev = (size_t *) malloc (MSBDF_MAX_ORD * sizeof (size_t));

  if (state->ordprev == 0)
    {
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for ordprev", GSL_ENOMEM);
    }

  state->ordprevbackup = (size_t *) malloc (MSBDF_MAX_ORD * sizeof (size_t));

  if (state->ordprevbackup == 0)
    {
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for ordprevbackup",
                      GSL_ENOMEM);
    }

  state->errlev = (double *) malloc (dim * sizeof (double));

  if (state->errlev == 0)
    {
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for errlev", GSL_ENOMEM);
    }

  state->abscor = gsl_vector_alloc (dim);

  if (state->abscor == 0)
    {
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for abscor", GSL_ENOMEM);
    }

  state->relcor = gsl_vector_alloc (dim);

  if (state->relcor == 0)
    {
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for relcor", GSL_ENOMEM);
    }

  state->svec = gsl_vector_alloc (dim);

  if (state->svec == 0)
    {
      gsl_vector_free (state->relcor);
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for svec", GSL_ENOMEM);
    }

  state->tempvec = gsl_vector_alloc (dim);

  if (state->tempvec == 0)
    {
      gsl_vector_free (state->svec);
      gsl_vector_free (state->relcor);
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for tempvec", GSL_ENOMEM);
    }

  state->dfdy = gsl_matrix_alloc (dim, dim);

  if (state->dfdy == 0)
    {
      gsl_vector_free (state->tempvec);
      gsl_vector_free (state->svec);
      gsl_vector_free (state->relcor);
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for dfdy", GSL_ENOMEM);
    }

  state->dfdt = (double *) malloc (dim * sizeof (double));

  if (state->dfdt == 0)
    {
      gsl_matrix_free (state->dfdy);
      gsl_vector_free (state->tempvec);
      gsl_vector_free (state->svec);
      gsl_vector_free (state->relcor);
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for dfdt", GSL_ENOMEM);
    }

  state->M = gsl_matrix_alloc (dim, dim);

  if (state->M == 0)
    {
      free (state->dfdt);
      gsl_matrix_free (state->dfdy);
      gsl_vector_free (state->tempvec);
      gsl_vector_free (state->svec);
      gsl_vector_free (state->relcor);
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for M", GSL_ENOMEM);
    }

  state->p = gsl_permutation_alloc (dim);

  if (state->p == 0)
    {
      gsl_matrix_free (state->M);
      free (state->dfdt);
      gsl_matrix_free (state->dfdy);
      gsl_vector_free (state->tempvec);
      gsl_vector_free (state->svec);
      gsl_vector_free (state->relcor);
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for p", GSL_ENOMEM);
    }

  state->rhs = gsl_vector_alloc (dim);

  if (state->rhs == 0)
    {
      gsl_permutation_free (state->p);
      gsl_matrix_free (state->M);
      free (state->dfdt);
      gsl_matrix_free (state->dfdy);
      gsl_vector_free (state->tempvec);
      gsl_vector_free (state->svec);
      gsl_vector_free (state->relcor);
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for rhs", GSL_ENOMEM);
    }

  state->abscorscaled = gsl_vector_alloc (dim);

  if (state->abscorscaled == 0)
    {
      gsl_vector_free (state->rhs); 
      gsl_permutation_free (state->p);
      gsl_matrix_free (state->M);
      free (state->dfdt);
      gsl_matrix_free (state->dfdy);
      gsl_vector_free (state->tempvec);
      gsl_vector_free (state->svec);
      gsl_vector_free (state->relcor);
      gsl_vector_free (state->abscor);
      free (state->errlev);
      free (state->ordprevbackup);
      free (state->ordprev);
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for abscorscaled", GSL_ENOMEM);
    }

  msbdf_reset ((void *) state, dim);

  state->driver = NULL;

  return state;
}

static int
msbdf_failurehandler (void *vstate, const size_t dim, const double t)
{
  /* Internal failure handler routine for msbdf. Adjusted strategy
     for GSL: Decrease order if this is the second time a failure
     has occurred at this order and point.
   */

  msbdf_state_t *state = (msbdf_state_t *) vstate;

  const size_t ord = state->ord;

  if (ord > 1 && (ord - state->ordprev[0] == 0) &&
      ord == state->failord && t == state->failt)
    {
      state->ord--;
    }

  /* Save information about failure */

  state->failord = ord;
  state->failt = t;
  state->ni++;

  /* Force reinitialization if failure took place at lowest
     order 
   */

  if (ord == 1)
    {
      msbdf_reset (vstate, dim);
    }

  return GSL_SUCCESS;
}

static int
msbdf_calccoeffs (const size_t ord, const size_t ordwait,
                  const double h, const double hprev[],
                  double l[],
                  double *errcoeff, double *ordm1coeff,
                  double *ordp1coeff, double *ordp2coeff, double *gamma)
{
  /* Calculates coefficients (l) of polynomial Lambda, error and
     auxiliary order change evaluation coefficients.
   */

  if (ord == 1)
    {
      l[0] = 1.0;
      l[1] = 1.0;
      *errcoeff = 0.5;
      *ordp1coeff = 2.0;

      {
        const double hsum = h + hprev[0];
        
        const double a5 = -1.5;
        const double a6 = -1.0 - h / hsum;
        const double c2 = 2.0 / (1.0 - a6 + a5);
        
        *ordp2coeff = fabs (c2 * (h / hsum) * 3.0 * a5);
      }
    }
  else
    {
      size_t i, j;
      double hsum = h;
      double coeff1 = -1.0;
      double x;

      /* Calculate the actual polynomial coefficients (l) */

      DBL_ZERO_MEMSET (l, MSBDF_MAX_ORD + 1);

      l[0] = 1.0;
      l[1] = 1.0;

      for (i = 2; i < ord; i++)
        {
          hsum += hprev[i - 2];
          coeff1 += -1.0 / i;

          for (j = i; j > 0; j--)
            {
              l[j] += h / hsum * l[j - 1];
            }
        }

      coeff1 += -1.0 / ord;

      x = -l[1] - coeff1;

      for (i = ord; i > 0; i--)
        {
          l[i] += l[i - 1] * x;
        }
      
#ifdef DEBUG
      {
        size_t di;
        
        printf ("-- calccoeffs l: ");
        for (di = 0; di < ord + 1; di++)
          {
            printf ("%.5e ", l[di]);
          }
        printf ("\n");
      }
#endif

      hsum += hprev[ord - 2];

      {
        const double coeff2 = -l[1] - h / hsum;
        const double a1 = 1.0 - coeff2 + coeff1;
        const double a2 = 1.0 + ord * a1;

        /* Calculate error coefficient */

        *errcoeff = fabs (a1 / (coeff1 * a2));

        /* Calculate auxiliary coefficients used in evaluation of change
           of order
        */

        if (ordwait < 2)
          {
            const double a3 = coeff1 + 1.0 / ord;
            const double a4 = coeff2 + h / hsum;
            const double c1 = a3 / (1.0 - a4 + a3);

            *ordm1coeff = fabs (c1 / (x / l[ord]));

            *ordp1coeff = fabs (a2 / (l[ord] * (h / hsum) / x));

            hsum += hprev[ord - 1];

            {
              const double a5 = coeff1 - 1.0 / (ord + 1.0);
              const double a6 = coeff2 - h / hsum;
              const double c2 = a2 / (1.0 - a6 + a5);

              *ordp2coeff = fabs (c2 * (h / hsum) * (ord + 2) * a5);
            }
          }
      }
    }

  *gamma = h / l[1];

#ifdef DEBUG
  printf ("-- calccoeffs ordm1coeff=%.5e ", *ordm1coeff);
  printf ("ordp1coeff=%.5e ", *ordp1coeff);
  printf ("ordp2coeff=%.5e ", *ordp2coeff);
  printf ("errcoeff=%.5e\n", *errcoeff);
#endif

  return GSL_SUCCESS;
}

static int
msbdf_update (void *vstate, const size_t dim, gsl_matrix * dfdy, double *dfdt,
              const double t, const double *y, const gsl_odeiv2_system * sys,
              gsl_matrix * M, gsl_permutation * p,
              const size_t iter, size_t * nJ, size_t * nM,
              const double tprev, const double failt,
              const double gamma, const double gammaprev, const double hratio)
{
  /* Evaluates Jacobian dfdy and updates iteration matrix M
     if criteria for update is met.
   */

  /* Jacobian is evaluated
     - at first step
     - if MSBDF_JAC_WAIT steps have been made without re-evaluation
     - in case of a convergence failure if
     --- change in gamma is small, or 
     --- convergence failure resulted in step size decrease
   */

  const double c = 0.2;
  const double gammarel = fabs (gamma / gammaprev - 1.0);

  if (*nJ == 0 || *nJ > MSBDF_JAC_WAIT ||
      (t == failt && (gammarel < c || hratio < 1.0)))
    {
#ifdef DEBUG
      printf ("-- evaluate jacobian\n");
#endif
      int s = GSL_ODEIV_JA_EVAL (sys, t, y, dfdy->data, dfdt);

      if (s == GSL_EBADFUNC)
        {
          return s;
        }

      if (s != GSL_SUCCESS)
        {
          msbdf_failurehandler (vstate, dim, t);
#ifdef DEBUG
          printf ("-- FAIL at jacobian function evaluation\n");
#endif
          return s;
        }

      /* Reset counter */

      *nJ = 0;
    }

  /* Iteration matrix M (and it's LU decomposition) is generated
     - at first step
     - if MSBDF_M_WAIT steps have been made without an update
     - if change in gamma is significant (e.g. change in step size)
     - if previous step was rejected
   */

  if (*nM == 0 || *nM > MSBDF_M_WAIT || gammarel >= c ||
      t == tprev || t == failt)
    {
#ifdef DEBUG
      printf ("-- update M, gamma=%.5e\n", gamma);
#endif
      size_t i;
      gsl_matrix_memcpy (M, dfdy);
      gsl_matrix_scale (M, -gamma);

      for (i = 0; i < dim; i++)
        {
          gsl_matrix_set (M, i, i, gsl_matrix_get (M, i, i) + 1.0);
        }

      {
        int signum;
        int s = gsl_linalg_LU_decomp (M, p, &signum);
        
        if (s != GSL_SUCCESS)
          {
            return GSL_FAILURE;
          }
      }

      /* Reset counter */

      *nM = 0;
    }

  return GSL_SUCCESS;
}

static int
msbdf_corrector (void *vstate, const gsl_odeiv2_system * sys,
                 const double t, const double h, const size_t dim,
                 const double z[], const double errlev[],
                 const double l[], const double errcoeff,
                 gsl_vector * abscor, gsl_vector * relcor,
                 double ytmp[], double ytmp2[],
                 gsl_matrix * dfdy, double dfdt[], gsl_matrix * M,
                 gsl_permutation * p, gsl_vector * rhs,
                 size_t * nJ, size_t * nM,
                 const double tprev, const double failt,
                 const double gamma, const double gammaprev,
                 const double hprev0)
{
  /* Calculates the correction step (abscor). Equation
     system M = I - gamma * dfdy = -G is solved by Newton iteration.
   */

  size_t mi, i;
  const size_t max_iter = 3;    /* Maximum number of iterations */
  double convrate = 1.0;        /* convergence rate */
  double stepnorm = 0.0;        /* norm of correction step */
  double stepnormprev = 0.0;    /* previous norm value */

  /* Evaluate at predicted values */

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h, z, ytmp);

    if (s == GSL_EBADFUNC)
      {
        return s;
      }

    if (s != GSL_SUCCESS)
      {
        msbdf_failurehandler (vstate, dim, t);

#ifdef DEBUG
        printf ("-- FAIL at user function evaluation\n");
#endif
        return s;
      }
  }

  /* Calculate correction step (abscor) */

  gsl_vector_set_zero (abscor);

  for (mi = 0; mi < max_iter; mi++)
    {
      const double safety = 0.3;
      const double safety2 = 0.1;

      /* Generate or update Jacobian and/or iteration matrix M if needed */

      if (mi == 0)
        {
          int s = msbdf_update (vstate, dim, dfdy, dfdt, t + h, z,
                                sys, M, p, mi,
                                nJ, nM, tprev, failt,
                                gamma, gammaprev,
                                h / hprev0);

          if (s != GSL_SUCCESS)
            {
              return s;
            }
        }

      /* Evaluate the right hand side (-G) */

      for (i = 0; i < dim; i++)
        {
          const double r = -1.0 * gsl_vector_get (abscor, i) -
            z[1 * dim + i] / l[1] + gamma * ytmp[i];

          gsl_vector_set (rhs, i, r);
        }

      /* Solve system of equations */

      {
        int s = gsl_linalg_LU_solve (M, p, rhs, relcor);
        
        if (s != GSL_SUCCESS)
          {
            msbdf_failurehandler (vstate, dim, t);
            
#ifdef DEBUG
            printf ("-- FAIL at LU_solve\n");
#endif
            return GSL_FAILURE;
          }
      }

#ifdef DEBUG
      {
        size_t di;
        printf ("-- dstep: ");
        for (di = 0; di < dim; di++)
          {
            printf ("%.5e ", gsl_vector_get (relcor, di));
          }
        printf ("\n");
      }
#endif

      /* Add iteration results */

      for (i = 0; i < dim; i++)
        {
          const double r =
            gsl_vector_get (abscor, i) + gsl_vector_get (relcor, i);

          gsl_vector_set (abscor, i, r);

          ytmp2[i] = z[i] + r;

          gsl_vector_set (relcor, i, gsl_vector_get (relcor, i) / errlev[i]);
        }

#ifdef DEBUG
      {
        size_t di;
        printf ("-- abscor: ");
        for (di = 0; di < dim; di++)
          {
            printf ("%.5e ", gsl_vector_get (abscor, di));
          }
        printf ("\n");
      }
#endif

      /* Convergence test. Norms used are root-mean-square norms. */

      stepnorm = gsl_blas_dnrm2 (relcor) / sqrt (dim);

      if (mi > 0)
        {
          convrate = GSL_MAX_DBL (safety * convrate, stepnorm / stepnormprev);
        }
      else
        {
          convrate = 1.0;
        }

      {
        const double convtest =
          GSL_MIN_DBL (convrate, 1.0) * stepnorm * errcoeff / safety2;

#ifdef DEBUG
        printf
          ("-- newt iter loop %d, errcoeff=%.5e, stepnorm =%.5e, convrate = %.5e, convtest = %.5e\n",
           (int) mi, errcoeff, stepnorm, convrate, convtest);
#endif
        if (convtest <= 1.0)
          {
            break;
          }
      }

      /* Check for divergence during iteration */

      {
        const double div_const = 2.0;

        if (mi > 1 && stepnorm > div_const * stepnormprev)
          {
            msbdf_failurehandler (vstate, dim, t);

#ifdef DEBUG
            printf ("-- FAIL, diverging Newton iteration\n");
#endif
            return GSL_FAILURE;
          }
      }

      /* Evaluate at new y */

      {
        int s = GSL_ODEIV_FN_EVAL (sys, t + h, ytmp2, ytmp);

        if (s == GSL_EBADFUNC)
          {
            return s;
          }

        if (s != GSL_SUCCESS)
          {
            msbdf_failurehandler (vstate, dim, t);

#ifdef DEBUG
            printf ("-- FAIL at user function evaluation\n");
#endif
            return s;
          }
      }

      stepnormprev = stepnorm;
    }

#ifdef DEBUG
  printf ("-- Newton iteration exit at mi=%d\n", (int) mi);
#endif

  /* Handle convergence failure */

  if (mi == max_iter)
    {
      msbdf_failurehandler (vstate, dim, t);

#ifdef DEBUG
      printf ("-- FAIL, max_iter reached\n");
#endif
      return GSL_FAILURE;
    }

  return GSL_SUCCESS;
}

static int
msbdf_eval_order (gsl_vector * abscorscaled, gsl_vector * tempvec,
                  gsl_vector * svec, const double errcoeff,
                  const size_t dim, const double errlev[],
                  const double ordm1coeff, const double ordp1coeff,
                  const double ordp1coeffprev, const double ordp2coeff,
                  const double hprev[],
                  const double h, const double z[],
                  size_t * ord)
{
  /* Evaluates and executes change in method order (current, current-1
     or current+1). Order which maximizes the step length is selected.
   */

  size_t i;

  /* step size estimates at current order, order-1 and order+1 */
  double ordest = 0.0;
  double ordm1est = 0.0;
  double ordp1est = 0.0;

  const double safety = 1e-6;
  const double bias = 6.0;
  const double bias2 = 10.0;
  const double min_incr = 1.5;

  /* Relative step length estimate for current order */

  ordest = 1.0 / (pow (bias * gsl_blas_dnrm2 (abscorscaled) / sqrt (dim)
                       * errcoeff, 1.0 / (*ord + 1)) + safety);

  /* Relative step length estimate for order ord - 1 */

  if (*ord > 1)
    {
      for (i = 0; i < dim; i++)
        {
          gsl_vector_set (tempvec, i, z[*ord * dim + i] / errlev[i]);
        }

      ordm1est = 1.0 / (pow (bias * gsl_blas_dnrm2 (tempvec) / sqrt (dim)
                             / ordm1coeff, 1.0 / (*ord)) + safety);
    }
  else
    {
      ordm1est = 0.0;
    }

  /* Relative step length estimate for order ord + 1 */

  if (*ord < MSBDF_MAX_ORD)
    {
      const double c = -ordp1coeff / ordp1coeffprev *
        pow (h / hprev[1], *ord + 1);

      for (i = 0; i < dim; i++)
        {
          gsl_vector_set (svec, i, gsl_vector_get (svec, i) * c +
                          gsl_vector_get (abscorscaled, i));
        }

      ordp1est = 1.0 / (pow (bias2 * gsl_blas_dnrm2 (svec) / sqrt (dim)
                             / ordp2coeff, 1.0 / (*ord + 2)) + safety);
    }
  else
    {
      ordp1est = 0.0;
    }

#ifdef DEBUG
  printf
    ("-- eval_order ord=%d, ordest=%.5e, ordm1est=%.5e, ordp1est=%.5e\n",
     (int) *ord, ordest, ordm1est, ordp1est);
#endif

  /* Choose order that maximises step size and increases step
     size markedly compared to current step 
   */

  if (ordm1est > ordest && ordm1est > ordp1est && ordm1est > min_incr)
    {
      *ord -= 1;
#ifdef DEBUG
      printf ("-- eval_order order DECREASED to %d\n", (int) *ord);
#endif
    }

  else if (ordp1est > ordest && ordp1est > ordm1est && ordp1est > min_incr)
    {
      *ord += 1;
#ifdef DEBUG
      printf ("-- eval_order order INCREASED to %d\n", (int) *ord);
#endif
    }

  return GSL_SUCCESS;
}

static int
msbdf_check_no_order_decrease (size_t const ordprev[])
{
  /* Checks if order has not been decreased according to order history
     array. Used in stability enhancement.
   */

  size_t i;

  for (i = 0; i < MSBDF_MAX_ORD - 1; i++)
    {
      if (ordprev[i + 1] > ordprev[i])
        {
          return 0;
        }
    }

  return 1;
}

static int
msbdf_check_step_size_decrease (double const hprev[])
{
  /* Checks if step size has decreased markedly according to
     step size history array. Used in stability enhancement.
   */

  size_t i;
  double max = fabs (hprev[0]);
  const double min = fabs (hprev[0]);
  const double decrease_limit = 0.5;

  for (i = 1; i < MSBDF_MAX_ORD; i++)
    {
      const double h = fabs (hprev[i]);

      if (h > min && h > max)
        {
          max = h;
        }
    }

  if (min / max < decrease_limit)
    {
      return 1;
    }

  return 0;
}

static int
msbdf_apply (void *vstate, size_t dim, double t, double h,
             double y[], double yerr[],
             const double dydt_in[], double dydt_out[],
             const gsl_odeiv2_system * sys)
{
  /* Carries out a step by BDF linear multistep methods. */

  msbdf_state_t *state = (msbdf_state_t *) vstate;

  double *const z = state->z;
  double *const zbackup = state->zbackup;
  double *const ytmp = state->ytmp;
  double *const ytmp2 = state->ytmp2;
  double *const l = state->l;
  double *const hprev = state->hprev;
  double *const hprevbackup = state->hprevbackup;
  size_t *const ordprev = state->ordprev;
  size_t *const ordprevbackup = state->ordprevbackup;
  double *const errlev = state->errlev;
  gsl_vector *const abscor = state->abscor;
  gsl_vector *const abscorscaled = state->abscorscaled;
  gsl_vector *const relcor = state->relcor;
  gsl_vector *const svec = state->svec;
  gsl_vector *const tempvec = state->tempvec;

  size_t ord = state->ord;      /* order for this step */
  double ordm1coeff = 0.0;
  double ordp1coeff = 0.0;
  double ordp2coeff = 0.0;
  double errcoeff = 0.0;        /* error coefficient */
  double gamma = 0.0;           /* gamma coefficient */

  const size_t max_failcount = 3;
  size_t i;

#ifdef DEBUG
  {
    size_t di;

    printf ("msbdf_apply: t=%.5e, ord=%d, h=%.5e, y:", t, (int) ord, h);

    for (di = 0; di < dim; di++)
      {
        printf ("%.5e ", y[di]);
      }
    printf ("\n");
  }
#endif

  /* Check if t is the same as on previous stepper call (or last
     failed call). This means that calculation of previous step failed
     or the step was rejected, and therefore previous state will be
     restored or the method will be reset.
   */

  if (state->ni > 0 && (t == state->tprev || t == state->failt))
    {
      if (state->ni == 1)
        {
          /* No step has been accepted yet, reset method */

          msbdf_reset (vstate, dim);
#ifdef DEBUG
          printf ("-- first step was REJECTED, msbdf_reset called\n");
#endif
        }
      else
        {
          /* A succesful step has been saved, restore previous state. */

          /* If previous step suggests order increase, but the step was
             rejected, then do not increase order.
           */

          if (ord > ordprev[0])
            {
              state->ord = ordprev[0];
              ord = state->ord;
            }

          /* Restore previous state */

          DBL_MEMCPY (z, zbackup, (MSBDF_MAX_ORD + 1) * dim);
          DBL_MEMCPY (hprev, hprevbackup, MSBDF_MAX_ORD);

          for (i = 0; i < MSBDF_MAX_ORD; i++)
            {
              ordprev[i] = ordprevbackup[i];
            }

          state->ordwait = state->ordwaitbackup;
          state->gammaprev = state->gammaprevbackup;

#ifdef DEBUG
          printf ("-- previous step was REJECTED, state restored\n");
#endif
        }

      /* If step is repeatedly rejected, then reset method */

      state->failcount++;

      if (state->failcount > max_failcount && state->ni > 1)
        {
          msbdf_reset (vstate, dim);
          ord = state->ord;

#ifdef DEBUG
          printf ("-- max_failcount reached, msbdf_reset called\n");
#endif
        }
    }
  else
    {
      /* The previous step was accepted. Backup current state. */

      DBL_MEMCPY (zbackup, z, (MSBDF_MAX_ORD + 1) * dim);
      DBL_MEMCPY (hprevbackup, hprev, MSBDF_MAX_ORD);

      for (i = 0; i < MSBDF_MAX_ORD; i++)
        {
          ordprevbackup[i] = ordprev[i];
        }

      state->ordwaitbackup = state->ordwait;
      state->gammaprevbackup = state->gammaprev;

      state->failcount = 0;

#ifdef DEBUG
      if (state->ni > 0)
        {
          printf ("-- previous step was ACCEPTED, state saved\n");
        }
#endif
    }

#ifdef DEBUG
  printf ("-- ord=%d, ni=%ld, ordwait=%d\n", (int) ord, state->ni,
          (int) state->ordwait);

  size_t di;
  printf ("-- ordprev: ");

  for (di = 0; di < MSBDF_MAX_ORD; di++)
    {
      printf ("%d ", (int) ordprev[di]);
    }

  printf ("\n");
#endif

  /* Get desired error levels via gsl_odeiv2_control object through driver
     object, which is a requirement for this stepper.
   */

  if (state->driver == NULL)
    {
      return GSL_EFAULT;
    }
  else
    {
      size_t i;

      for (i = 0; i < dim; i++)
        {
          if (dydt_in != NULL)
            {
              gsl_odeiv2_control_errlevel (state->driver->c, y[i],
                                           dydt_in[i], h, i, &errlev[i]);
            }
          else
            {
              gsl_odeiv2_control_errlevel (state->driver->c, y[i],
                                           0.0, h, i, &errlev[i]);
            }
        }
    }

#ifdef DEBUG
  {
    size_t di;
    printf ("-- errlev: ");
    for (di = 0; di < dim; di++)
      {
        printf ("%.5e ", errlev[di]);
      }
    printf ("\n");
  }
#endif

  /* On first call initialize Nordsieck matrix */

  if (state->ni == 0)
    {
      size_t i;

      DBL_ZERO_MEMSET (z, (MSBDF_MAX_ORD + 1) * dim);

      if (dydt_in != NULL)
        {
          DBL_MEMCPY (ytmp, dydt_in, dim);
        }
      else
        {
          int s = GSL_ODEIV_FN_EVAL (sys, t, y, ytmp);

          if (s != GSL_SUCCESS)
            {
              return s;
            }
        }

      DBL_MEMCPY (&z[0 * dim], y, dim);
      DBL_MEMCPY (&z[1 * dim], ytmp, dim);

      for (i = 0; i < dim; i++)
        {
          z[1 * dim + i] *= h;
        }
    }

  /* Stability enhancement heuristic for msbdf: If order > 1 and order
     has not been changed, check for decrease in step size, that is
     not accompanied by a decrease in method order. This condition may
     be indication of BDF method stability problems, a change in ODE
     system, or convergence problems in Newton iteration. In all
     cases, the strategy is to decrease method order.
   */

#ifdef DEBUG
  printf ("-- check_no_order_decrease %d, check_step_size_decrease %d\n",
          msbdf_check_no_order_decrease (ordprev),
          msbdf_check_step_size_decrease (hprev));
#endif

  if (ord > 1 &&
      ord - ordprev[0] == 0 &&
      msbdf_check_no_order_decrease (ordprev) &&
      msbdf_check_step_size_decrease (hprev))
    {
      state->ord--;
      state->ordwait = ord + 2;
      ord = state->ord;

#ifdef DEBUG
      printf ("-- stability enhancement decreased order to %d\n", (int) ord);
#endif
    }

  /* Sanity check */

  { 
    const int deltaord = ord - ordprev[0];

  if (deltaord > 1 || deltaord < -1)
    {
      printf ("-- order change %d\n", deltaord);
      GSL_ERROR ("msbdf_apply too large order change", GSL_ESANITY);
    }

  /* Modify Nordsieck matrix if order or step length has been changed */

  /* If order increased by 1, adjust Nordsieck matrix */

  if (deltaord == 1)
    {
      if (ord > 2)
        {
          size_t i, j;
          double hsum = h;
          double coeff1 = -1.0;
          double coeff2 = 1.0;
          double hrelprev = 1.0;
          double hrelprod = 1.0;
          double hrel = 0.0;

          /* Calculate coefficients used in adjustment to l */

          DBL_ZERO_MEMSET (l, MSBDF_MAX_ORD + 1);

          l[2] = 1.0;

          for (i = 1; i < ord - 1; i++)
            {
              hsum += hprev[i];
              hrel = hsum / h;
              hrelprod *= hrel;
              coeff1 -= 1.0 / (i + 1);
              coeff2 += 1.0 / hrel;

              for (j = i + 2; j > 1; j--)
                {
                  l[j] *= hrelprev;
                  l[j] += l[j - 1];
                }

              hrelprev = hrel;
            }

          /* Scale Nordsieck matrix */

          {
            const double c = (-coeff1 - coeff2) / hrelprod;

            for (i = 0; i < dim; i++)
              {
                z[ord * dim + i] = c * gsl_vector_get (abscor, i);
              }
          }
          for (i = 2; i < ord; i++)
            for (j = 0; j < dim; j++)
              {
                z[i * dim + j] += l[i] * z[ord * dim + j];
              }
        }
      else
        {
          /* zero new vector for order incease from 1 to 2 */

          DBL_ZERO_MEMSET (&z[ord * dim], dim);
        }

#ifdef DEBUG
      printf ("-- order increase detected, Nordsieck modified\n");
#endif
    }

  /* If order decreased by 1, adjust Nordsieck matrix */

  if (deltaord == -1)
    {
      size_t i, j;
      double hsum = 0.0;

      /* Calculate coefficients used in adjustment to l */

      DBL_ZERO_MEMSET (l, MSBDF_MAX_ORD + 1);

      l[2] = 1.0;

      for (i = 1; i < ord; i++)
        {
          hsum += hprev[i - 1];

          for (j = i + 2; j > 1; j--)
            {
              l[j] *= hsum / h;
              l[j] += l[j - 1];
            }
        }

      /* Scale Nordsieck matrix */

      for (i = 2; i < ord + 1; i++)
        for (j = 0; j < dim; j++)
          {
            z[i * dim + j] += -l[i] * z[(ord + 1) * dim + j];
          }

#ifdef DEBUG
      printf ("-- order decrease detected, Nordsieck modified\n");
#endif
    }

  /* Scale Nordsieck vectors if step size has been changed */

  if (state->ni > 0 && h != hprev[0])
    {
      size_t i, j;
      const double hrel = h / hprev[0];
      double coeff = hrel;

      for (i = 1; i < ord + 1; i++)
        {
          for (j = 0; j < dim; j++)
            {
              z[i * dim + j] *= coeff;
            }

          coeff *= hrel;
        }

#ifdef DEBUG
      printf ("-- h != hprev, Nordsieck modified\n");
#endif
    }

  /* Calculate polynomial coefficients (l), error coefficient and
     auxiliary coefficients
   */

  msbdf_calccoeffs (ord, state->ordwait, h, hprev, l, &errcoeff,
                    &ordm1coeff, &ordp1coeff, &ordp2coeff, &gamma);

  /* Carry out the prediction step */

  {
    size_t i, j, k;

    for (i = 1; i < ord + 1; i++)
      for (j = ord; j > i - 1; j--)
        for (k = 0; k < dim; k++)
          {
            z[(j - 1) * dim + k] += z[j * dim + k];
          }

#ifdef DEBUG
    {
      size_t di;
      printf ("-- predicted y: ");
      for (di = 0; di < dim; di++)
        {
          printf ("%.5e ", z[di]);
        }
      printf ("\n");
    }
#endif
  }

  /* Calculate correction step to abscor */
  {
    int s;
    s = msbdf_corrector (vstate, sys, t, h, dim, z, errlev, l, errcoeff,
                         abscor, relcor, ytmp, ytmp2,
                         state->dfdy, state->dfdt, state->M,
                         state->p, state->rhs,
                         &(state->nJ), &(state->nM),
                         state->tprev, state->failt, gamma,
                         state->gammaprev, hprev[0]);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  {
    /* Add accepted final correction step to Nordsieck matrix */

    size_t i, j;

    for (i = 0; i < ord + 1; i++)
      for (j = 0; j < dim; j++)
        {
          z[i * dim + j] += l[i] * gsl_vector_get (abscor, j);
        }

#ifdef DEBUG
    {
      size_t di;
      printf ("---- l: ");
      for (di = 0; di < ord + 1; di++)
        {
          printf ("%.5e ", l[di]);
        }
      printf ("\n");

      printf ("-- corrected y: ");
      for (di = 0; di < dim; di++)
        {
          printf ("%.5e ", z[di]);
        }
      printf ("\n");
    }
#endif

    /* Derivatives at output */

    if (dydt_out != NULL)
      {
        int s = GSL_ODEIV_FN_EVAL (sys, t + h, z, dydt_out);

        if (s == GSL_EBADFUNC)
          {
            return s;
          }

        if (s != GSL_SUCCESS)
          {
            msbdf_failurehandler (vstate, dim, t);

#ifdef DEBUG
            printf ("-- FAIL at user function evaluation\n");
#endif
            return s;
          }
      }

    /* Calculate error estimate */

    for (i = 0; i < dim; i++)
      {
        yerr[i] = fabs (gsl_vector_get (abscor, i)) * errcoeff;
      }

#ifdef DEBUG
    {
      size_t di;
      printf ("-- yerr: ");
      for (di = 0; di < dim; di++)
        {
          printf ("%.5e ", yerr[di]);
        }
      printf ("\n");
    }
#endif

    /* Save y values */

    for (i = 0; i < dim; i++)
      {
        y[i] = z[0 * dim + i];
      }
  }

  /* Scale abscor with errlev for later use in norm calculations
     in order evaluation in msbdf_eval_order
   */
  {
    size_t i;

    for (i = 0; i < dim; i++)
      {
        gsl_vector_set (abscorscaled, i, 
                        gsl_vector_get (abscor, i) / errlev[i]);
      }
  }

  /* Save items needed for evaluation of order increase on next
     call, if needed
   */

  if (state->ordwait == 1 && ord < MSBDF_MAX_ORD)
    {
      size_t i;

      state->ordp1coeffprev = ordp1coeff;

      for (i = 0; i < dim; i++)
        {
          gsl_vector_set (svec, i, gsl_vector_get (abscorscaled, i));
        }
    }

  /* Consider and execute order change for next step if order is unchanged. */

  if (state->ordwait == 0) 
  { 
    if (ord - ordprev[0] == 0)
      {
        msbdf_eval_order (abscorscaled, tempvec, svec, errcoeff, dim, errlev,
                          ordm1coeff, ordp1coeff,
                          state->ordp1coeffprev, ordp2coeff,
                          hprev, h, z, &(state->ord));

        state->ordwait = state->ord + 2;
      }
    else
      {
        /* Postpone order evaluation if order has been modified elsewhere */
    
        state->ordwait = 2;
      }
  }
  
  /* Save information about current step in state and update counters */
  {
    size_t i;
    
    for (i = MSBDF_MAX_ORD - 1; i > 0; i--)
      {
        hprev[i] = hprev[i - 1];
        ordprev[i] = ordprev[i - 1];
      }
  }
  
  hprev[0] = h;
  ordprev[0] = ord;
  
#ifdef DEBUG
  {
    size_t di;
    printf ("-- hprev: ");
    for (di = 0; di < MSBDF_MAX_ORD; di++)
      {
        printf ("%.5e ", hprev[di]);
      }
    printf ("\n");
  }
#endif
  
  state->tprev = t;
  state->ordwait--;
  state->ni++;
  state->gammaprev = gamma;
  
  state->nJ++;
  state->nM++;
  
#ifdef DEBUG
  printf ("-- nJ=%d, nM=%d\n", (int) state->nJ, (int) state->nM);
#endif
  }

  return GSL_SUCCESS;
}

static int
msbdf_set_driver (void *vstate, const gsl_odeiv2_driver * d)
{
  msbdf_state_t *state = (msbdf_state_t *) vstate;

  state->driver = d;

  return GSL_SUCCESS;
}

static int
msbdf_reset (void *vstate, size_t dim)
{
  msbdf_state_t *state = (msbdf_state_t *) vstate;
  size_t i;

  state->ni = 0;
  state->ord = 1;
  state->ordwait = 2;
  state->ordwaitbackup = 2;
  state->failord = 0;
  state->failt = GSL_NAN;
  state->tprev = GSL_NAN;  
  state->gammaprev = 1.0;
  state->gammaprevbackup = 1.0;
  state->nJ = 0;
  state->nM = 0;
  state->failcount = 0;
  state->ordp1coeffprev = 0.0;

  DBL_ZERO_MEMSET (state->hprev, MSBDF_MAX_ORD);
  DBL_ZERO_MEMSET (state->hprevbackup, MSBDF_MAX_ORD);
  DBL_ZERO_MEMSET (state->z, (MSBDF_MAX_ORD + 1) * dim);
  DBL_ZERO_MEMSET (state->zbackup, (MSBDF_MAX_ORD + 1) * dim);

  for (i = 0; i < MSBDF_MAX_ORD; i++)
    {
      state->ordprev[i] = 1;
      state->ordprevbackup[i] = 1;
    }

#ifdef DEBUG
  printf ("-- msbdf_reset called\n");
#endif

  return GSL_SUCCESS;
}

static unsigned int
msbdf_order (void *vstate)
{
  msbdf_state_t *state = (msbdf_state_t *) vstate;

  return state->ord;
}

static void
msbdf_free (void *vstate)
{
  msbdf_state_t *state = (msbdf_state_t *) vstate;

  gsl_vector_free (state->rhs);
  gsl_permutation_free (state->p);
  gsl_matrix_free (state->M);
  free (state->dfdt);
  gsl_matrix_free (state->dfdy);
  gsl_vector_free (state->tempvec);
  gsl_vector_free (state->svec);
  gsl_vector_free (state->relcor);
  gsl_vector_free (state->abscor);
  gsl_vector_free (state->abscorscaled);
  free (state->errlev);
  free (state->ordprevbackup);
  free (state->ordprev);
  free (state->hprevbackup);
  free (state->hprev);
  free (state->l);
  free (state->ytmp2);
  free (state->ytmp);
  free (state->zbackup);
  free (state->z);
  free (state);
}

static const gsl_odeiv2_step_type msbdf_type = {
  "msbdf",                      /* name */
  1,                            /* can use dydt_in? */
  1,                            /* gives exact dydt_out? */
  &msbdf_alloc,
  &msbdf_apply,
  &msbdf_set_driver,
  &msbdf_reset,
  &msbdf_order,
  &msbdf_free
};

const gsl_odeiv2_step_type *gsl_odeiv2_step_msbdf = &msbdf_type;
