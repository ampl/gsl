/* ode-initval2/msadams.c
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

/* A variable-coefficient linear multistep Adams method in Nordsieck
   form. This stepper uses explicit Adams-Bashforth (predictor) and
   implicit Adams-Moulton (corrector) methods in P(EC)^m functional
   iteration mode. Method order varies dynamically between 1 and 12.

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

#include "odeiv_util.h"

/* Maximum order of Adams methods */
#define MSADAMS_MAX_ORD 12

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
  double *pc;                   /* product term coefficients */
  double *l;                    /* polynomial coefficients */
  double *hprev;                /* previous step sizes */
  double *hprevbackup;          /* backup of hprev */
  double *errlev;               /* desired error level of y */
  gsl_vector *abscor;           /* absolute y values for correction */
  gsl_vector *relcor;           /* relative y values for correction */
  gsl_vector *svec;             /* saved abscor & work area */
  gsl_vector *tempvec;          /* work area */
  const gsl_odeiv2_driver *driver;      /* pointer to gsl_odeiv2_driver object */

  long int ni;                  /* stepper call counter */
  size_t ord;                   /* current order of method */
  size_t ordprev;               /* order of previous call */
  size_t ordprevbackup;         /* backup of ordprev */
  double tprev;                 /* t point of previous call */
  size_t ordwait;               /* counter for order change */
  size_t ordwaitbackup;         /* backup of ordwait */
  size_t failord;               /* order of convergence failure */
  double failt;                 /* t point of convergence failure */
  double ordm1coeff;            /* saved order-1 coefficiet */
  double ordp1coeffprev;        /* saved order+1 coefficient */
  size_t failcount;             /* counter for rejected steps */
}
msadams_state_t;

/* Introduce msadams_reset for use in msadams_alloc and _apply */

static int msadams_reset (void *, size_t);

static void *
msadams_alloc (size_t dim)
{
  msadams_state_t *state =
    (msadams_state_t *) malloc (sizeof (msadams_state_t));

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for msadams_state",
                      GSL_ENOMEM);
    }

  state->z =
    (double *) malloc ((MSADAMS_MAX_ORD + 1) * dim * sizeof (double));

  if (state->z == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for z", GSL_ENOMEM);
    }

  state->zbackup =
    (double *) malloc ((MSADAMS_MAX_ORD + 1) * dim * sizeof (double));

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

  state->pc = (double *) malloc ((MSADAMS_MAX_ORD + 1) * sizeof (double));

  if (state->pc == 0)
    {
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for pc", GSL_ENOMEM);
    }

  state->l = (double *) malloc ((MSADAMS_MAX_ORD + 1) * sizeof (double));

  if (state->l == 0)
    {
      free (state->pc);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for l", GSL_ENOMEM);
    }

  state->hprev = (double *) malloc (MSADAMS_MAX_ORD * sizeof (double));

  if (state->hprev == 0)
    {
      free (state->l);
      free (state->pc);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for hprev", GSL_ENOMEM);
    }

  state->hprevbackup = (double *) malloc (MSADAMS_MAX_ORD * sizeof (double));

  if (state->hprevbackup == 0)
    {
      free (state->hprev);
      free (state->l);
      free (state->pc);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for hprevbackup", GSL_ENOMEM);
    }

  state->errlev = (double *) malloc (dim * sizeof (double));

  if (state->errlev == 0)
    {
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->pc);
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
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->pc);
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
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->pc);
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
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->pc);
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
      free (state->hprevbackup);
      free (state->hprev);
      free (state->l);
      free (state->pc);
      free (state->ytmp2);
      free (state->ytmp);
      free (state->zbackup);
      free (state->z);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for tempvec", GSL_ENOMEM);
    }

  msadams_reset ((void *) state, dim);

  state->driver = NULL;

  return state;
}

static int
msadams_failurehandler (void *vstate, const size_t dim, const double t)
{
  /* Internal failure handler routine for msadams. Adjusted strategy
     for GSL: Decrease order if this is the second time a failure
     has occurred at this order and point.
   */

  msadams_state_t *state = (msadams_state_t *) vstate;

  const size_t ord = state->ord;

  if (ord > 1 && (ord - state->ordprev == 0) &&
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
      msadams_reset (vstate, dim);
    }

  return GSL_SUCCESS;
}

static int
msadams_calccoeffs (const size_t ord, const size_t ordwait,
                    const double h, const double hprev[],
                    double pc[], double l[],
                    double *errcoeff, double *ordm1coeff,
                    double *ordp1coeff, double *ordp2coeff)
{

  /* Calculates coefficients (l) of polynomial Lambda, error and
     auxiliary order change evaluation coefficients.
   */

  if (ord == 1)
    {
      l[0] = 1.0;
      l[1] = 1.0;
      *errcoeff = 0.5;
      *ordp1coeff = 1.0;
      *ordp2coeff = 12.0;
    }
  else
    {
      size_t i, j;
      double hsum = h;
      double st1 = 0.0;         /* sum term coefficients */
      double st2 = 0.0;

      /* Calculate coefficients (pc) of product terms */

      DBL_ZERO_MEMSET (pc, MSADAMS_MAX_ORD + 1);

      pc[0] = 1.0;

      for (i = 1; i < ord; i++)
        {
          /* Calculate auxiliary coefficient used in evaluation of
             change of order
           */

          if (i == ord - 1 && ordwait == 1)
            {
              int s = 1;

              *ordm1coeff = 0.0;

              for (j = 0; j < ord - 1; j++)
                {
                  *ordm1coeff += s * pc[j] / (j + 2);
                  s = -s;
                }

              *ordm1coeff = pc[ord - 2] / (ord * (*ordm1coeff));
            }

          for (j = i; j > 0; j--)
            {
              pc[j] += pc[j - 1] * h / hsum;
            }

          hsum += hprev[i - 1];
        }

      /* Calculate sum term 1 for error estimation */

      {
        int s = 1;

        for (i = 0; i < ord; i++)
          {
            st1 += s * pc[i] / (i + 1.0);
            s = -s;
          }
      }

      /* Calculate sum term 2 for error estimation */

      {
        int s = 1;

        for (i = 0; i < ord; i++)
          {
            st2 += s * pc[i] / (i + 2.0);
            s = -s;
          }
      }

      /* Calculate the actual polynomial coefficients (l) */

      DBL_ZERO_MEMSET (l, MSADAMS_MAX_ORD + 1);

      l[0] = 1.0;

      for (i = 1; i < ord + 1; i++)
        {
          l[i] = pc[i - 1] / (i * st1);
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

        printf ("-- calccoeffs pc: ");
        for (di = 0; di < ord; di++)
          {
            printf ("%.5e ", pc[di]);
          }
        printf ("\n");

        printf ("-- calccoeffs st1=%.5e, st2=%.5e\n", st1, st2);
      }
#endif

      /* Calculate error coefficient */

      *errcoeff = (h / hsum) * (st2 / st1);

      /* Calculate auxiliary coefficients used in evaluation of change
         of order
       */

      if (ordwait < 2)
        {
          int s = 1;

          *ordp1coeff = hsum / (h * l[ord]);

          *ordp2coeff = 0.0;

          for (i = ord; i > 0; i--)
            {
              pc[i] += pc[i - 1] * (h / hsum);
            }

          for (i = 0; i < ord + 1; i++)
            {
              *ordp2coeff += s * pc[i] / (i + 2);
              s = -s;
            }

          *ordp2coeff = (ord + 1) * st1 / (*ordp2coeff);
        }
    }

#ifdef DEBUG
  printf ("-- calccoeffs ordm1coeff=%.5e ", *ordm1coeff);
  printf ("ordp1coeff=%.5e ", *ordp1coeff);
  printf ("ordp2coeff=%.5e ", *ordp2coeff);
  printf ("errcoeff=%.5e\n", *errcoeff);
#endif

  return GSL_SUCCESS;
}

static int
msadams_corrector (void *vstate, const gsl_odeiv2_system * sys,
                   const double t, const double h, const size_t dim,
                   const double z[], const double errlev[],
                   const double l[], const double errcoeff,
                   gsl_vector * abscor, gsl_vector * relcor,
                   double ytmp[], double ytmp2[])
{
  /* Calculates the correction step (abscor). Non-linear equation
     system is solved by functional iteration.
   */

  size_t mi, i;
  const size_t max_iter = 3;    /* Maximum number of iterations */
  double convrate = 1.0;        /* convergence rate */
  double stepnorm = 0.0;        /* norm of correction step */
  double stepnormprev = 0.0;

  /* Evaluate at predicted values */

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h, z, ytmp);

    if (s == GSL_EBADFUNC)
      {
        return s;
      }

    if (s != GSL_SUCCESS)
      {
        msadams_failurehandler (vstate, dim, t);

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

      /* Calculate new y values to ytmp2 */

      for (i = 0; i < dim; i++)
        {
          ytmp[i] *= h;
          ytmp[i] -= z[1 * dim + i];
          ytmp[i] /= l[1];
          ytmp2[i] = z[0 * dim + i] + ytmp[i];
        }

#ifdef DEBUG
      {
        size_t di;
        printf ("-- dstep: ");
        for (di = 0; di < dim; di++)
          {
            printf ("%.5e ", ytmp[di]);
          }
        printf ("\n");
      }
#endif
      /* Convergence test. Norms used are root-mean-square norms. */

      for (i = 0; i < dim; i++)
        {
          gsl_vector_set (relcor, i,
                          (ytmp[i] - gsl_vector_get (abscor, i)) / errlev[i]);
          gsl_vector_set (abscor, i, ytmp[i]);
        }

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
          ("-- func iter loop %d, errcoeff=%.5e, stepnorm =%.5e, convrate = %.5e, convtest = %.5e\n",
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
            msadams_failurehandler (vstate, dim, t);
            
#ifdef DEBUG
            printf ("-- FAIL, diverging functional iteration\n");
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
            msadams_failurehandler (vstate, dim, t);

#ifdef DEBUG
            printf ("-- FAIL at user function evaluation\n");
#endif
            return s;
          }
      }

      stepnormprev = stepnorm;
    }

#ifdef DEBUG
  printf ("-- functional iteration exit at mi=%d\n", (int) mi);
#endif

  /* Handle convergence failure */

  if (mi == max_iter)
    {
      msadams_failurehandler (vstate, dim, t);

#ifdef DEBUG
      printf ("-- FAIL, max_iter reached\n");
#endif
      return GSL_FAILURE;
    }

  return GSL_SUCCESS;
}

static int
msadams_eval_order (gsl_vector * abscor, gsl_vector * tempvec,
                    gsl_vector * svec, const double errcoeff,
                    const size_t dim, const double errlev[],
                    const double ordm1coeff, const double ordp1coeff,
                    const double ordp1coeffprev, const double ordp2coeff,
                    const double hprev[],
                    const double h, const double z[],
                    size_t * ord, size_t * ordwait)
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

  /* Relative step length estimate for current order */

  ordest = 1.0 / (pow (bias * gsl_blas_dnrm2 (abscor) / sqrt (dim)
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

  if (*ord < MSADAMS_MAX_ORD)
    {
      const double c = -ordp1coeff / ordp1coeffprev *
        pow (h / hprev[1], *ord + 1);

      for (i = 0; i < dim; i++)
        {
          gsl_vector_set (svec, i, gsl_vector_get (svec, i) * c +
                          gsl_vector_get (abscor, i));
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

  {
    const double min_incr = 1.5;

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
  }

  *ordwait = *ord + 2;


  return GSL_SUCCESS;
}


static int
msadams_apply (void *vstate, size_t dim, double t, double h,
               double y[], double yerr[],
               const double dydt_in[], double dydt_out[],
               const gsl_odeiv2_system * sys)
{
  /* Conducts a step by Adams linear multistep methods. */

  msadams_state_t *state = (msadams_state_t *) vstate;

  double *const z = state->z;
  double *const zbackup = state->zbackup;
  double *const ytmp = state->ytmp;
  double *const ytmp2 = state->ytmp2;
  double *const pc = state->pc;
  double *const l = state->l;
  double *const hprev = state->hprev;
  double *const hprevbackup = state->hprevbackup;
  double *const errlev = state->errlev;
  gsl_vector *const abscor = state->abscor;
  gsl_vector *const relcor = state->relcor;
  gsl_vector *const svec = state->svec;
  gsl_vector *const tempvec = state->tempvec;

  size_t ord = state->ord;
  double ordm1coeff = 0.0;
  double ordp1coeff = 0.0;
  double ordp2coeff = 0.0;
  double errcoeff = 0.0;        /* error coefficient */

  int deltaord;

#ifdef DEBUG
  {
    size_t di;

    printf ("msadams_apply: t=%.5e, ord=%d, h=%.5e, y:", t, (int) ord, h);

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

          msadams_reset (vstate, dim);
#ifdef DEBUG
          printf ("-- first step was REJECTED, msadams_reset called\n");
#endif
        }
      else
        {
          /* A succesful step has been saved, restore previous state. */

          /* If previous step suggests order increase, but the step was
             rejected, then do not increase order.
           */

          if (ord > state->ordprev)
            {
              state->ord = state->ordprev;
              ord = state->ord;
            }

          /* Restore previous state */

          DBL_MEMCPY (z, zbackup, (MSADAMS_MAX_ORD + 1) * dim);
          DBL_MEMCPY (hprev, hprevbackup, MSADAMS_MAX_ORD);
          state->ordprev = state->ordprevbackup;
          state->ordwait = state->ordwaitbackup;

#ifdef DEBUG
          printf ("-- previous step was REJECTED, state restored\n");
#endif
        }

      state->failcount++;

      {
        const size_t max_failcount = 3;
        
        /* If step is repeatedly rejected, then reset method */

        if (state->failcount > max_failcount && state->ni > 1)
          {
            msadams_reset (vstate, dim);
            ord = state->ord;
            
#ifdef DEBUG
            printf ("-- max_failcount reached, msadams_reset called\n");
#endif
          }

        /* Attempt to decrease order twice is indicative of system stiffness.
           Reset method to enable fast recovery. */

        else if ((int)state->ordprev - (int)ord >= 2)
          {
            msadams_reset (vstate, dim);
            ord = state->ord;
          
#ifdef DEBUG
            printf ("-- order decreased by two, msadams_reset called\n");
#endif
          } 
      }
    }
  else
    {
      /* The previous step was accepted. Backup current state. */

      DBL_MEMCPY (zbackup, z, (MSADAMS_MAX_ORD + 1) * dim);
      DBL_MEMCPY (hprevbackup, hprev, MSADAMS_MAX_ORD);
      state->ordprevbackup = state->ordprev;
      state->ordwaitbackup = state->ordwait;

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
  printf ("-- ordprev: %d\n", (int) state->ordprev);
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

      DBL_ZERO_MEMSET (z, (MSADAMS_MAX_ORD + 1) * dim);

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

  /* Sanity check  */

  deltaord = ord - state->ordprev;

  if (deltaord > 1 || deltaord < -1)
    {
      printf ("-- order change %d\n", deltaord);
      GSL_ERROR ("msadams_apply too large order change", GSL_ESANITY);
    }

  /* Modify Nordsieck matrix if order or step length has been changed */

  /* If order increased by 1, initialize new Nordsieck vector */

  if (deltaord == 1)
    {
      DBL_ZERO_MEMSET (&z[ord * dim], dim);

#ifdef DEBUG
      printf ("-- order increase detected, Nordsieck modified\n");
#endif
    }

  /* If order decreased by 1, adjust Nordsieck matrix */

  if (deltaord == -1)
    {
      double hsum = 0.0;
      size_t i, j;

      /* Calculate coefficients used in adjustment to l */

      DBL_ZERO_MEMSET (l, MSADAMS_MAX_ORD + 1);

      l[1] = 1.0;

      for (i = 1; i < ord; i++)
        {
          hsum += hprev[i - 1];

          for (j = i + 1; j > 0; j--)
            {
              l[j] *= hsum / hprev[0];
              l[j] += l[j - 1];
            }
        }

      for (i = 1; i < ord; i++)
        {
          l[i + 1] = (ord + 1) * l[i] / (i + 1);
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

  msadams_calccoeffs (ord, state->ordwait, h, hprev, pc, l, &errcoeff,
                      &ordm1coeff, &ordp1coeff, &ordp2coeff);

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
    s = msadams_corrector (vstate, sys, t, h, dim, z, errlev, l, errcoeff,
                           abscor, relcor, ytmp, ytmp2);
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
            msadams_failurehandler (vstate, dim, t);

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

  /* Scale abscor with errlev for later use in norm calculations */
  {
    size_t i;

    for (i = 0; i < dim; i++)
      {
        gsl_vector_set (abscor, i, gsl_vector_get (abscor, i) / errlev[i]);
      }
  }

  /* Save items needed for evaluation of order increase on next
     call, if needed
   */

  if (state->ordwait == 1 && ord < MSADAMS_MAX_ORD)
    {
      size_t i;

      state->ordp1coeffprev = ordp1coeff;
      state->ordm1coeff = ordm1coeff;

      for (i = 0; i < dim; i++)
        {
          gsl_vector_set (svec, i, gsl_vector_get (abscor, i));
        }
    }

  /* Consider and execute order change for next step */

  if (state->ordwait == 0)
    {
      msadams_eval_order (abscor, tempvec, svec, errcoeff, dim, errlev,
                          state->ordm1coeff, ordp1coeff,
                          state->ordp1coeffprev, ordp2coeff,
                          hprev, h, z, &(state->ord), &(state->ordwait));
    }

  /* Save information about current step in state and update counters */
  {
    size_t i;

    state->ordprev = ord;

    for (i = MSADAMS_MAX_ORD - 1; i > 0; i--)
      {
        hprev[i] = hprev[i - 1];
      }
    hprev[0] = h;

#ifdef DEBUG
    {
      size_t di;
      printf ("-- hprev: ");
      for (di = 0; di < MSADAMS_MAX_ORD; di++)
        {
          printf ("%.5e ", hprev[di]);
        }
      printf ("\n");
    }
#endif

    state->tprev = t;
    state->ordwait--;
    state->ni++;
  }

  return GSL_SUCCESS;
}

static int
msadams_set_driver (void *vstate, const gsl_odeiv2_driver * d)
{
  msadams_state_t *state = (msadams_state_t *) vstate;

  state->driver = d;

  return GSL_SUCCESS;
}

static int
msadams_reset (void *vstate, size_t dim)
{
  msadams_state_t *state = (msadams_state_t *) vstate;

  state->ni = 0;
  state->ord = 1;
  state->ordprev = 1;
  state->ordprevbackup = 1;
  state->ordwait = 2;
  state->ordwaitbackup = 2;
  state->failord = 0;
  state->failt = GSL_NAN;
  state->failcount = 0;

  DBL_ZERO_MEMSET (state->hprev, MSADAMS_MAX_ORD);
  DBL_ZERO_MEMSET (state->z, (MSADAMS_MAX_ORD + 1) * dim);

#ifdef DEBUG
  printf ("-- msadams_reset called\n");
#endif

  return GSL_SUCCESS;
}

static unsigned int
msadams_order (void *vstate)
{
  msadams_state_t *state = (msadams_state_t *) vstate;

  return state->ord;
}

static void
msadams_free (void *vstate)
{
  msadams_state_t *state = (msadams_state_t *) vstate;

  gsl_vector_free (state->tempvec);
  gsl_vector_free (state->svec);
  gsl_vector_free (state->relcor);
  gsl_vector_free (state->abscor);
  free (state->errlev);
  free (state->hprevbackup);
  free (state->hprev);
  free (state->l);
  free (state->pc);
  free (state->ytmp2);
  free (state->ytmp);
  free (state->zbackup);
  free (state->z);
  free (state);
}

static const gsl_odeiv2_step_type msadams_type = {
  "msadams",                    /* name */
  1,                            /* can use dydt_in? */
  1,                            /* gives exact dydt_out? */
  &msadams_alloc,
  &msadams_apply,
  &msadams_set_driver,
  &msadams_reset,
  &msadams_order,
  &msadams_free
};

const gsl_odeiv2_step_type *gsl_odeiv2_step_msadams = &msadams_type;
