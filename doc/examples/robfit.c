#include <stdio.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>

int
dofit(const gsl_multifit_robust_type *T,
      const gsl_matrix *X, const gsl_vector *y,
      gsl_vector *c, gsl_matrix *cov)
{
  int s;
  gsl_multifit_robust_workspace * work 
    = gsl_multifit_robust_alloc (T, X->size1, X->size2);

  s = gsl_multifit_robust (X, y, c, cov, work);
  gsl_multifit_robust_free (work);

  return s;
}

int
main (int argc, char **argv)
{
  int i;
  size_t n;
  const size_t p = 2; /* linear fit */
  gsl_matrix *X, *cov;
  gsl_vector *x, *y, *c, *c_ols;
  const double a = 1.45; /* slope */
  const double b = 3.88; /* intercept */
  gsl_rng *r;

  if (argc != 2)
    {
      fprintf (stderr,"usage: robfit n\n");
      exit (-1);
    }

  n = atoi (argv[1]);

  X = gsl_matrix_alloc (n, p);
  x = gsl_vector_alloc (n);
  y = gsl_vector_alloc (n);

  c = gsl_vector_alloc (p);
  c_ols = gsl_vector_alloc (p);
  cov = gsl_matrix_alloc (p, p);

  r = gsl_rng_alloc(gsl_rng_default);

  /* generate linear dataset */
  for (i = 0; i < n - 3; i++)
    {
      double dx = 10.0 / (n - 1.0);
      double ei = gsl_rng_uniform(r);
      double xi = -5.0 + i * dx;
      double yi = a * xi + b;

      gsl_vector_set (x, i, xi);
      gsl_vector_set (y, i, yi + ei);
    }

  /* add a few outliers */
  gsl_vector_set(x, n - 3, 4.7);
  gsl_vector_set(y, n - 3, -8.3);

  gsl_vector_set(x, n - 2, 3.5);
  gsl_vector_set(y, n - 2, -6.7);

  gsl_vector_set(x, n - 1, 4.1);
  gsl_vector_set(y, n - 1, -6.0);

  /* construct design matrix X for linear fit */
  for (i = 0; i < n; ++i)
    {
      double xi = gsl_vector_get(x, i);

      gsl_matrix_set (X, i, 0, 1.0);
      gsl_matrix_set (X, i, 1, xi);
    }

  /* perform robust and OLS fit */
  dofit(gsl_multifit_robust_ols, X, y, c_ols, cov);
  dofit(gsl_multifit_robust_bisquare, X, y, c, cov);

  /* output data and model */
  for (i = 0; i < n; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(y, i);
      gsl_vector_view v = gsl_matrix_row(X, i);
      double y_ols, y_rob, y_err;

      gsl_multifit_robust_est(&v.vector, c, cov, &y_rob, &y_err);
      gsl_multifit_robust_est(&v.vector, c_ols, cov, &y_ols, &y_err);

      printf("%g %g %g %g\n", xi, yi, y_rob, y_ols);
    }

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  {
    printf ("# best fit: Y = %g + %g X\n", 
            C(0), C(1));

    printf ("# covariance matrix:\n");
    printf ("# [ %+.5e, %+.5e\n",
               COV(0,0), COV(0,1));
    printf ("#   %+.5e, %+.5e\n", 
               COV(1,0), COV(1,1));
  }

  gsl_matrix_free (X);
  gsl_vector_free (x);
  gsl_vector_free (y);
  gsl_vector_free (c);
  gsl_vector_free (c_ols);
  gsl_matrix_free (cov);
  gsl_rng_free(r);

  return 0;
}
