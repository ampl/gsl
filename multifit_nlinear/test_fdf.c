/* multifit_nlinear/test_fdf.c
 * 
 * Copyright (C) 2007, 2013, 2014, 2015, 2016 Brian Gough, Patrick Alken
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

typedef struct
{
  const char *name;
  double *x0;       /* initial parameters (size p) */
  double *weights;  /* data weights */
  double *sigma;
  double *epsrel;   /* relative tolerance for solution checking */
  void (*checksol) (const double x[], const double sumsq,
                    const double epsrel, const char *sname,
                    const char *pname);
  gsl_multifit_nlinear_fdf *fdf;
} test_fdf_problem;

#include "test_bard.c"
#include "test_beale.c"
#include "test_biggs.c"
#include "test_box.c"
#include "test_boxbod.c"
#include "test_brown1.c"
#include "test_brown2.c"
#include "test_brown3.c"
#include "test_eckerle.c"
#include "test_enso.c"
#include "test_exp1.c"
#include "test_gaussian.c"
#include "test_hahn1.c"
#include "test_helical.c"
#include "test_jennrich.c"
#include "test_kirby2.c"
#include "test_kowalik.c"
#include "test_lin1.c"
#include "test_lin2.c"
#include "test_lin3.c"
#include "test_meyer.c"
#include "test_meyerscal.c"
#include "test_osborne.c"
#include "test_penalty1.c"
#include "test_penalty2.c"
#include "test_powell1.c"
#include "test_powell2.c"
#include "test_powell3.c"
#include "test_rat42.c"
#include "test_rat43.c"
#include "test_rosenbrock.c"
#include "test_rosenbrocke.c"
#include "test_roth.c"
#include "test_thurber.c"
#include "test_vardim.c"
#include "test_watson.c"
#include "test_wood.c"

#include "test_wnlin.c"

static void test_fdf(const gsl_multifit_nlinear_type * T,
                     const gsl_multifit_nlinear_parameters * params,
                     const double xtol, const double gtol,
                     const double ftol,
                     const double epsrel,
                     test_fdf_problem *problem);
static void test_fdf_checksol(const char *sname, const char *pname,
                              const double epsrel,
                              gsl_multifit_nlinear_workspace *s,
                              test_fdf_problem *problem);

/*
 * FIXME: some test problems are disabled since they fail on certain
 * solvers. Known failures are:
 *
 * Method     test-problem
 * ======     ============
 * dogleg     thurbera
 * dogleg     rat43a
 * all        boxboda
 */

static test_fdf_problem *test_problems[] = {
  /*
   * These test problems are taken from
   *
   * H. B. Nielsen, UCTP test problems for unconstrained optimization,
   * IMM Department of Mathematical Modeling, Tech. Report
   * IMM-REP-2000-17, 2000.
   */
  &lin1_problem,       /* 1 */
  &lin2_problem,       /* 2 */
  &lin3_problem,       /* 3 */
  &rosenbrock_problem, /* 4 */
  &helical_problem,    /* 5 */
  &powell1_problem,    /* 6 */
  &roth_problem,       /* 7 */
  &bard_problem,       /* 8 */
  &kowalik_problem,    /* 9 */
  &meyer_problem,      /* 10 */
  &watson_problem,     /* 11 */
  &box_problem,        /* 12 */
  &jennrich_problem,   /* 13 */
  &brown1_problem,     /* 14 */
  &brown2_problem,     /* 16 */
  &osborne_problem,    /* 17 */
  &exp1_problem,       /* 18 */
  &meyerscal_problem,  /* 20 */

  &powell2_problem,

  /*
   * These tests are from
   *
   * J. J. More, B. S. Garbow and K. E. Hillstrom, Testing
   * Unconstrained Optimization Software, ACM Trans. Math. Soft.
   * Vol 7, No 1, 1981.
   *
   * Many of these overlap with the Nielsen tests
   */
  &rosenbrock_problem,   /* 1 */
  &roth_problem,         /* 2 */
  &powell3_problem,      /* 3 */
  &brown3_problem,       /* 4 */
  &beale_problem,        /* 5 */
  &jennrich_problem,     /* 6 */
  &helical_problem,      /* 7 */
  &bard_problem,         /* 8 */
  &gaussian_problem,     /* 9 */
  &meyer_problem,        /* 10 */
  &box_problem,          /* 12 */
  &powell1_problem,      /* 13 */
  &wood_problem,         /* 14 */
  &kowalik_problem,      /* 15 */
  &brown1_problem,       /* 16 */
  &osborne_problem,      /* 17 */
  &biggs_problem,        /* 18 */
  &watson_problem,       /* 20 */
  &rosenbrocke_problem,  /* 21 */
  &penalty1_problem,     /* 23 */
  &penalty2_problem,     /* 24 */
  &vardim_problem,       /* 25 */
  &brown2_problem,       /* 27 */
  &lin1_problem,         /* 32 */
  &lin2_problem,         /* 33 */
  &lin3_problem,         /* 34 */

  /* NIST test cases */
  &kirby2a_problem,
  &kirby2b_problem,
  &hahn1a_problem,
  &hahn1b_problem,
  &ensoa_problem,
  &ensob_problem,
  /*&thurbera_problem,*/
  &thurberb_problem,
  /*&boxboda_problem,*/
  &boxbodb_problem,
  &rat42a_problem,
  &rat42b_problem,
  &eckerlea_problem,
  &eckerleb_problem,
  /*&rat43a_problem,*/
  &rat43b_problem,

  /* weighted test cases */
  &wnlin_problem1,
  &wnlin_problem2,

  NULL
};

static void
test_fdf_main(const gsl_multifit_nlinear_parameters * params)
{
  const double xtol = pow(GSL_DBL_EPSILON, 0.9);
  const double gtol = pow(GSL_DBL_EPSILON, 0.9);
  const double ftol = 0.0;
  size_t i;

  for (i = 0; test_problems[i] != NULL; ++i)
    {
      test_fdf_problem *problem = test_problems[i];
      double epsrel = *(problem->epsrel);
      gsl_multifit_nlinear_fdf fdf;

      test_fdf(gsl_multifit_nlinear_trust, params, xtol, gtol, ftol,
               epsrel, problem);

      /* test finite difference Jacobian
       * XXX: watson problem doesn't work with forward differences */
      if (problem != &watson_problem)
        {
          fdf.df = problem->fdf->df;
          problem->fdf->df = NULL;

          test_fdf(gsl_multifit_nlinear_trust, params, xtol, gtol, ftol,
                   1.0e3 * epsrel, problem);

          problem->fdf->df = fdf.df;
        }

#if 0 /*XXX: box3d test fails on MacOS here */
      if (params->trs == gsl_multifit_nlinear_trs_lmaccel && problem->fdf->fvv != NULL)
        {
          /* test finite difference second directional derivative */
          fdf.fvv = problem->fdf->fvv;
          problem->fdf->fvv = NULL;

          test_fdf(gsl_multifit_nlinear_trust, params, xtol, gtol, ftol,
                   epsrel / params->h_fvv, problem);

          problem->fdf->fvv = fdf.fvv;
        }
#endif
    }
}

/*
test_fdf()
  Test a weighted nonlinear least squares problem

Inputs: T        - solver to use
        params   - solver parameters
        xtol     - tolerance in x
        gtol     - tolerance in gradient
        ftol     - tolerance in residual vector
        epsrel   - relative error tolerance in solution
        problem  - contains the nonlinear problem and solution point
*/

static void
test_fdf(const gsl_multifit_nlinear_type * T,
         const gsl_multifit_nlinear_parameters * params,
         const double xtol, const double gtol, const double ftol,
         const double epsrel,
         test_fdf_problem *problem)
{
  gsl_multifit_nlinear_fdf *fdf = problem->fdf;
  const size_t n = fdf->n;
  const size_t p = fdf->p;
  const size_t max_iter = 2500;
  gsl_vector *x0 = gsl_vector_alloc(p);
  gsl_vector_view x0v = gsl_vector_view_array(problem->x0, p);
  gsl_multifit_nlinear_workspace *w =
    gsl_multifit_nlinear_alloc (T, params, n, p);
  const char *pname = problem->name;
  char buf[1024];
  char sname[2048];
  int status, info;

  sprintf(buf, "%s/%s/scale=%s/solver=%s",
    gsl_multifit_nlinear_name(w),
    params->trs->name,
    params->scale->name,
    params->solver->name);

  if (problem->fdf->df == NULL)
    {
      if (params->fdtype == GSL_MULTIFIT_NLINEAR_FWDIFF)
        strcat(buf, "/fdjac,forward");
      else
        strcat(buf, "/fdjac,center");
    }

  if (problem->fdf->fvv == NULL)
    {
      strcat(buf, "/fdfvv");
    }

  strcpy(sname, buf);

  gsl_vector_memcpy(x0, &x0v.vector);

  if (problem->weights != NULL)
    {
      gsl_vector_const_view wv = gsl_vector_const_view_array(problem->weights, n);
      gsl_multifit_nlinear_winit(x0, &wv.vector, fdf, w);
    }
  else
    gsl_multifit_nlinear_init(x0, fdf, w);

  status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                                       NULL, NULL, &info, w);
  gsl_test(status, "%s/%s did not converge, status=%s",
           sname, pname, gsl_strerror(status));

  /* check solution */
  test_fdf_checksol(sname, pname, epsrel, w, problem);

  if (problem->weights == NULL)
    {
      /* test again with weighting matrix W = I */
      gsl_vector *wv = gsl_vector_alloc(n);

      sprintf(sname, "%s/weighted", buf);

      gsl_vector_memcpy(x0, &x0v.vector);

      gsl_vector_set_all(wv, 1.0);
      gsl_multifit_nlinear_winit(x0, wv, fdf, w);
  
      status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                                           NULL, NULL, &info, w);
      gsl_test(status, "%s/%s did not converge, status=%s",
               sname, pname, gsl_strerror(status));

      test_fdf_checksol(sname, pname, epsrel, w, problem);

      gsl_vector_free(wv);
    }

  gsl_multifit_nlinear_free(w);
  gsl_vector_free(x0);
}

static void
test_fdf_checksol(const char *sname, const char *pname,
                  const double epsrel,
                  gsl_multifit_nlinear_workspace *w,
                  test_fdf_problem *problem)
{
  gsl_multifit_nlinear_fdf *fdf = problem->fdf;
  const double *sigma = problem->sigma;
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double sumsq;

  /* check solution vector x and sumsq = ||f||^2 */
  gsl_blas_ddot(f, f, &sumsq);
  (problem->checksol)(x->data, sumsq, epsrel, sname, pname);

  /* check variances */
  if (sigma)
    {
      const size_t n = fdf->n;
      const size_t p = fdf->p;
      size_t i;
      gsl_matrix * covar = gsl_matrix_alloc (p, p);
      gsl_matrix *J = gsl_multifit_nlinear_jac (w);

      gsl_multifit_nlinear_covar (J, 0.0, covar);

      for (i = 0; i < p; i++) 
        {
          double ei = sqrt(sumsq/(n-p))*sqrt(gsl_matrix_get(covar,i,i));
          gsl_test_rel (ei, sigma[i], epsrel, 
                        "%s/%s, sigma(%d)", sname, pname, i) ;
        }

      gsl_matrix_free (covar);
    }
}
