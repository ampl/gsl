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
  const double b = 2.0 * M_PI;
  const size_t spline_order = 6;     /* spline order */
  const size_t ncontrol = 15;        /* number of control points */
  const double sigma = 0.2;          /* noise */
  gsl_bspline_workspace *w = gsl_bspline_alloc_ncontrol(spline_order, ncontrol);
  gsl_bspline_workspace *wper = gsl_bspline_alloc_ncontrol(spline_order, ncontrol);
  gsl_vector *c = gsl_vector_alloc(ncontrol);    /* non-periodic coefficients */
  gsl_vector *cper = gsl_vector_alloc(ncontrol); /* periodic coefficients */
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *wts = gsl_vector_alloc(n);
  size_t i;
  gsl_rng *r;
  double chisq;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);

  /* this is the data to be fitted */
  for (i = 0; i < n; ++i)
    {
      double xi = (b - a) / (n - 1.0) * i + a;
      double yi = sin(xi) - cos(2.0 * xi);
      double dyi = gsl_ran_gaussian(r, sigma);

      yi += dyi;

      gsl_vector_set(x, i, xi);
      gsl_vector_set(y, i, yi);
      gsl_vector_set(wts, i, 1.0 / (sigma * sigma));

      printf("%f %f\n", xi, yi);
    }

  printf("\n\n");

  /* use uniform non-periodic knots on [a, b] */
  gsl_bspline_init_uniform(a, b, w);

  /* solve least squares problem for non-periodic spline */
  gsl_bspline_wlssolve(x, y, wts, c, &chisq, w);

  /* use periodic knots on [a, b] */
  gsl_bspline_init_periodic(a, b, wper);

  /* solve least squares problem for periodic spline */
  gsl_bspline_pwlssolve(x, y, wts, cper, &chisq, wper);

  /* output the spline curves */
  {
    double xi;

    for (xi = a; xi <= b; xi += 0.01)
      {
        double result, result_per;
        gsl_bspline_calc(xi, c, &result, w);
        gsl_bspline_calc(xi, cper, &result_per, wper);
        printf("%f %f %f\n", xi, result, result_per);
      }
  }

  fprintf(stderr, "=== Non-periodic spline endpoint derivatives ===\n");

  for (i = 0; i < spline_order; ++i)
    {
      double result0, result1;

      gsl_bspline_calc_deriv(a, c, i, &result0, w);
      gsl_bspline_calc_deriv(b, c, i, &result1, w);

      fprintf(stderr, "deriv %zu: [%14.6e, %14.6e]\n", i, result0, result1);
    }

  fprintf(stderr, "=== Periodic spline endpoint derivatives ===\n");

  for (i = 0; i < spline_order; ++i)
    {
      double result0, result1;

      gsl_bspline_calc_deriv(a, cper, i, &result0, wper);
      gsl_bspline_calc_deriv(b, cper, i, &result1, wper);

      fprintf(stderr, "deriv %zu: [%14.6e, %14.6e]\n", i, result0, result1);
    }

  gsl_vector_free(c);
  gsl_vector_free(cper);
  gsl_bspline_free(w);
  gsl_bspline_free(wper);
  gsl_rng_free(r);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(wts);

  return 0;
}
