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
  const size_t K = 7;                                  /* window size */
  const double f = 5.0;                                /* frequency of square wave in Hz */
  gsl_filter_median_workspace *median_p = gsl_filter_median_alloc(K);
  gsl_filter_rmedian_workspace *rmedian_p = gsl_filter_rmedian_alloc(K);
  gsl_vector *t = gsl_vector_alloc(N);                 /* time */
  gsl_vector *x = gsl_vector_alloc(N);                 /* input vector */
  gsl_vector *y_median = gsl_vector_alloc(N);          /* median filtered output */
  gsl_vector *y_rmedian = gsl_vector_alloc(N);         /* recursive median filtered output */
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  size_t i;

  /* generate input signal */
  for (i = 0; i < N; ++i)
    {
      double ti = (double) i / (N - 1.0);
      double tmp = sin(2.0 * M_PI * f * ti);
      double xi = (tmp >= 0.0) ? 1.0 : -1.0;
      double ei = gsl_ran_gaussian(r, 0.1);

      gsl_vector_set(t, i, ti);
      gsl_vector_set(x, i, xi + ei);
    }

  gsl_filter_median(GSL_FILTER_END_PADVALUE, x, y_median, median_p);
  gsl_filter_rmedian(GSL_FILTER_END_PADVALUE, x, y_rmedian, rmedian_p);

  /* print results */
  for (i = 0; i < N; ++i)
    {
      double ti = gsl_vector_get(t, i);
      double xi = gsl_vector_get(x, i);
      double medi = gsl_vector_get(y_median, i);
      double rmedi = gsl_vector_get(y_rmedian, i);

      printf("%f %f %f %f\n",
             ti,
             xi,
             medi,
             rmedi);
    }

  gsl_vector_free(t);
  gsl_vector_free(x);
  gsl_vector_free(y_median);
  gsl_vector_free(y_rmedian);
  gsl_rng_free(r);
  gsl_filter_median_free(median_p);

  return 0;
}
