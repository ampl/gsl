#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

int
main(void)
{
  const size_t N = 500;                      /* length of time series */
  const size_t K = 11;                       /* window size */
  gsl_movstat_workspace * w = gsl_movstat_alloc(K);
  gsl_vector *x = gsl_vector_alloc(N);
  gsl_vector *xmean = gsl_vector_alloc(N);
  gsl_vector *xmin = gsl_vector_alloc(N);
  gsl_vector *xmax = gsl_vector_alloc(N);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double xi = cos(4.0 * M_PI * i / (double) N);
      double ei = gsl_ran_gaussian(r, 0.1);

      gsl_vector_set(x, i, xi + ei);
    }

  /* compute moving statistics */
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE, x, xmean, w);
  gsl_movstat_minmax(GSL_MOVSTAT_END_PADVALUE, x, xmin, xmax, w);

  /* print results */
  for (i = 0; i < N; ++i)
    {
      printf("%zu %f %f %f %f\n",
             i,
             gsl_vector_get(x, i),
             gsl_vector_get(xmean, i),
             gsl_vector_get(xmin, i),
             gsl_vector_get(xmax, i));
    }

  gsl_vector_free(x);
  gsl_vector_free(xmean);
  gsl_rng_free(r);
  gsl_movstat_free(w);

  return 0;
}
