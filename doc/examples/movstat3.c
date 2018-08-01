#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

double
func(const size_t n, double x[], void * params)
{
  const double alpha = *(double *) params;

  gsl_sort(x, 1, n);

  return gsl_stats_trmean_from_sorted_data(alpha, x, 1, n);
}

int
main(void)
{
  const size_t N = 1000;                /* length of time series */
  const size_t K = 11;                  /* window size */
  double alpha = 0.1;                   /* trimmed mean parameter */
  gsl_vector *x = gsl_vector_alloc(N);  /* input vector */
  gsl_vector *y = gsl_vector_alloc(N);  /* filtered output vector for alpha1 */
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_movstat_workspace *w = gsl_movstat_alloc(K);
  gsl_movstat_function F;
  size_t i;
  double sum = 0.0;

  /* generate input signal */
  for (i = 0; i < N; ++i)
    {
      double ui = gsl_ran_gaussian(r, 1.0);
      double outlier = (gsl_rng_uniform(r) < 0.01) ? 10.0*GSL_SIGN(ui) : 0.0;
      sum += ui;
      gsl_vector_set(x, i, sum + outlier);
    }

  /* apply moving window function */
  F.function = func;
  F.params = &alpha;
  gsl_movstat_apply(GSL_MOVSTAT_END_PADVALUE, &F, x, y, w);

  /* print results */
  for (i = 0; i < N; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(y, i);

      printf("%f %f\n", xi, yi);
    }

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_rng_free(r);
  gsl_movstat_free(w);

  return 0;
}
