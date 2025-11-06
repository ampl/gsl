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
  const size_t n = 500;              /* number of data points to fit */
  const double a = 0.0;              /* data interval [a,b] */
  const double b = 15.0;
  const size_t max_spline_order = 5; /* maximum spline order */
  const size_t nbreak = 10;          /* number of breakpoints */
  const double sigma = 0.2;          /* noise */
  gsl_bspline_workspace **work = malloc(max_spline_order * sizeof(gsl_bspline_workspace *));
  gsl_vector **c = malloc(max_spline_order * sizeof(gsl_vector *));
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *w = gsl_vector_alloc(n);
  size_t i;
  gsl_rng *r;
  double chisq;

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

  printf("\n\n");

  for (i = 0; i < max_spline_order; ++i)
    {
      /* allocate workspace for this spline order */
      work[i] = gsl_bspline_alloc(i + 1, nbreak);
      c[i] = gsl_vector_alloc(gsl_bspline_ncontrol(work[i]));

      /* use uniform breakpoints on [a, b] */
      gsl_bspline_init_uniform(a, b, work[i]);

      /* solve least squares problem */
      gsl_bspline_wlssolve(x, y, w, c[i], &chisq, work[i]);
    }

  /* output the spline curves */
  {
    double xi;

    for (xi = a; xi <= b; xi += 0.1)
      {
        printf("%f ", xi);
        for (i = 0; i < max_spline_order; ++i)
          {
            double result;
            gsl_bspline_calc(xi, c[i], &result, work[i]);
            printf("%f ", result);
          }
        printf("\n");
      }
  }

  for (i = 0; i < max_spline_order; ++i)
    {
      gsl_vector_free(c[i]);
      gsl_bspline_free(work[i]);
    }

  gsl_rng_free(r);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(w);
  free(work);
  free(c);

  return 0;
}
