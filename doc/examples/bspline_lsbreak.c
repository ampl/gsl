#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

int
main (void)
{
  const size_t n = 500;     /* number of data points to fit */
  const size_t k = 4;       /* spline order */
  const double a = 0.0;     /* data interval [a,b] */
  const double b = 15.0;
  const double sigma = 0.2; /* noise */
  gsl_bspline_workspace *work1 = gsl_bspline_alloc(k, 40);   /* 40 breakpoints */
  gsl_bspline_workspace *work2 = gsl_bspline_alloc(k, 10);   /* 10 breakpoints */
  const size_t p1 = gsl_bspline_ncontrol(work1);             /* number of control points */
  const size_t p2 = gsl_bspline_ncontrol(work2);             /* number of control points */
  const size_t dof1 = n - p1;                                /* degrees of freedom */
  const size_t dof2 = n - p2;                                /* degrees of freedom */
  gsl_vector *c1 = gsl_vector_alloc(p1);
  gsl_vector *c2 = gsl_vector_alloc(p2);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *w = gsl_vector_alloc(n);
  size_t i;
  gsl_rng *r;
  double chisq1, chisq2;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);

  /* this is the data to be fitted */
  for (i = 0; i < n; ++i)
    {
      double xi = (b - a) / (n - 1.0) * i + a;
      double yi = cos(xi) * exp(-0.1 * xi);
      double dyi = gsl_ran_gaussian(r, sigma);

      yi += dyi;

      gsl_vector_set(x, i, xi);
      gsl_vector_set(y, i, yi);
      gsl_vector_set(w, i, 1.0 / (sigma * sigma));

      printf("%f %f\n", xi, yi);
    }

  /* use uniform breakpoints on [a, b] */
  gsl_bspline_init_uniform(a, b, work1);
  gsl_bspline_init_uniform(a, b, work2);

  /* solve least squares problem */
  gsl_bspline_wlssolve(x, y, w, c1, &chisq1, work1);
  gsl_bspline_wlssolve(x, y, w, c2, &chisq2, work2);

  fprintf(stderr, "40 breakpoints: chisq/dof = %e\n", chisq1 / dof1);
  fprintf(stderr, "10 breakpoints: chisq/dof = %e\n", chisq2 / dof2);

  printf("\n\n");

  /* output the spline curves */
  {
    double xi;

    for (xi = a; xi <= b; xi += 0.1)
      {
        double result1, result2;

        gsl_bspline_calc(xi, c1, &result1, work1);
        gsl_bspline_calc(xi, c2, &result2, work2);
        printf("%f %f %f\n", xi, result1, result2);
      }
  }

  gsl_rng_free(r);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c1);
  gsl_vector_free(c2);
  gsl_bspline_free(work1);
  gsl_bspline_free(work2);

  return 0;
}
