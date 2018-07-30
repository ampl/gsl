#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

int
main(void)
{
  const size_t N = 1000;                       /* length of time series */
  const size_t K = 61;                         /* window size */
  const double alpha = 3.0;                    /* Gaussian kernel has +/- 3 standard deviations */
  gsl_vector *x = gsl_vector_alloc(N);         /* input vector */
  gsl_vector *y = gsl_vector_alloc(N);         /* filtered output vector */
  gsl_vector *dy = gsl_vector_alloc(N);        /* first derivative filtered vector */
  gsl_vector *d2y = gsl_vector_alloc(N);       /* second derivative filtered vector */
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_filter_gaussian_workspace *gauss_p = gsl_filter_gaussian_alloc(K);
  size_t i;

  /* generate input signal */
  for (i = 0; i < N; ++i)
    {
      double xi = (i > N / 2) ? 0.5 : 0.0;
      double ei = gsl_ran_gaussian(r, 0.1);

      gsl_vector_set(x, i, xi + ei);
    }

  /* apply filters */
  gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha, 0, x, y, gauss_p);
  gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha, 1, x, dy, gauss_p);
  gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha, 2, x, d2y, gauss_p);

  /* print results */
  for (i = 0; i < N; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(y, i);
      double dyi = gsl_vector_get(dy, i);
      double d2yi = gsl_vector_get(d2y, i);
      double dxi;

      /* compute finite difference of x vector */
      if (i == 0)
        dxi = gsl_vector_get(x, i + 1) - xi;
      else if (i == N - 1)
        dxi = gsl_vector_get(x, i) - gsl_vector_get(x, i - 1);
      else
        dxi = 0.5 * (gsl_vector_get(x, i + 1) - gsl_vector_get(x, i - 1));

      printf("%.12e %.12e %.12e %.12e %.12e\n",
             xi,
             yi,
             dxi,
             dyi,
             d2yi);
    }

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(dy);
  gsl_vector_free(d2y);
  gsl_rng_free(r);
  gsl_filter_gaussian_free(gauss_p);

  return 0;
}
