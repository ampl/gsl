#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "expfit.c"

/* number of data points to fit */
#define N 40

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver *s;
  int status, info;
  size_t i;
  const size_t n = N;
  const size_t p = 3;

  gsl_matrix *J = gsl_matrix_alloc(n, p);
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double y[N], weights[N];
  struct data d = { n, y };
  gsl_multifit_function_fdf f;
  double x_init[3] = { 1.0, 0.0, 0.0 };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  gsl_vector_view w = gsl_vector_view_array(weights, n);
  const gsl_rng_type * type;
  gsl_rng * r;
  gsl_vector *res_f;
  double chi, chi0;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &expb_f;
  f.df = &expb_df;   /* set to NULL for finite-difference Jacobian */
  f.n = n;
  f.p = p;
  f.params = &d;

  /* This is the data to be fitted */

  for (i = 0; i < n; i++)
    {
      double t = i;
      double yi = 1.0 + 5 * exp (-0.1 * t);
      double si = 0.1 * yi;
      double dy = gsl_ran_gaussian(r, si);

      weights[i] = 1.0 / (si * si);
      y[i] = yi + dy;
      printf ("data: %zu %g %g\n", i, y[i], si);
    };

  s = gsl_multifit_fdfsolver_alloc (T, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_fdfsolver_wset (s, &f, &x.vector, &w.vector);

  /* compute initial residual norm */
  res_f = gsl_multifit_fdfsolver_residual(s);
  chi0 = gsl_blas_dnrm2(res_f);

  /* solve the system with a maximum of 20 iterations */
  status = gsl_multifit_fdfsolver_driver(s, 20, xtol, gtol, ftol, &info);

  gsl_multifit_fdfsolver_jac(s, J);
  gsl_multifit_covar (J, 0.0, covar);

  /* compute final residual norm */
  chi = gsl_blas_dnrm2(res_f);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  fprintf(stderr, "summary from method '%s'\n",
          gsl_multifit_fdfsolver_name(s));
  fprintf(stderr, "number of iterations: %zu\n",
          gsl_multifit_fdfsolver_niter(s));
  fprintf(stderr, "function evaluations: %zu\n", f.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", f.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",
          (info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial |f(x)| = %g\n", chi0);
  fprintf(stderr, "final   |f(x)| = %g\n", chi);

  { 
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

    fprintf(stderr, "chisq/dof = %g\n",  pow(chi, 2.0) / dof);

    fprintf (stderr, "A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    fprintf (stderr, "lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    fprintf (stderr, "b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
  }

  fprintf (stderr, "status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_matrix_free (J);
  gsl_rng_free (r);
  return 0;
}
