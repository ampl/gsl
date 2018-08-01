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
  const size_t N = 1000;                               /* length of time series */
  const size_t K = 25;                                 /* window size */
  const double t = 4.0;                                /* number of scale factors for outlier detection */
  gsl_vector *x = gsl_vector_alloc(N);                 /* input vector */
  gsl_vector *y = gsl_vector_alloc(N);                 /* output (filtered) vector */
  gsl_vector *xmedian = gsl_vector_alloc(N);           /* window medians */
  gsl_vector *xsigma = gsl_vector_alloc(N);            /* window scale estimates */
  gsl_vector_int *ioutlier = gsl_vector_int_alloc(N);  /* outlier detected? */
  gsl_filter_impulse_workspace * w = gsl_filter_impulse_alloc(K);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  size_t noutlier;
  size_t i;

  /* generate input signal */
  for (i = 0; i < N; ++i)
    {
      double xi = 10.0 * sin(2.0 * M_PI * i / (double) N);
      double ei = gsl_ran_gaussian(r, 2.0);
      double u = gsl_rng_uniform(r);
      double outlier = (u < 0.01) ? 15.0*GSL_SIGN(ei) : 0.0;

      gsl_vector_set(x, i, xi + ei + outlier);
    }

  /* apply impulse detection filter */
  gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_QN, t, x, y,
                     xmedian, xsigma, &noutlier, ioutlier, w);

  /* print results */
  for (i = 0; i < N; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(y, i);
      double xmedi = gsl_vector_get(xmedian, i);
      double xsigmai = gsl_vector_get(xsigma, i);
      int outlier = gsl_vector_int_get(ioutlier, i);

      printf("%zu %f %f %f %f %d\n",
             i,
             xi,
             yi,
             xmedi + t * xsigmai,
             xmedi - t * xsigmai,
             outlier);
    }

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(xmedian);
  gsl_vector_free(xsigma);
  gsl_vector_int_free(ioutlier);
  gsl_filter_impulse_free(w);
  gsl_rng_free(r);

  return 0;
}
