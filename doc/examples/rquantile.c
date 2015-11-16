#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>

int
main(void)
{
  const size_t N = 10000;
  double *data = malloc(N * sizeof(double));
  gsl_rstat_quantile_workspace *work_25 = gsl_rstat_quantile_alloc(0.25);
  gsl_rstat_quantile_workspace *work_50 = gsl_rstat_quantile_alloc(0.5);
  gsl_rstat_quantile_workspace *work_75 = gsl_rstat_quantile_alloc(0.75);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  double exact_p25, exact_p50, exact_p75;
  double val_p25, val_p50, val_p75;
  size_t i;

  /* add data to quantile accumulators; also store data for exact
   * comparisons */
  for (i = 0; i < N; ++i)
    {
      data[i] = gsl_ran_rayleigh(r, 1.0);
      gsl_rstat_quantile_add(data[i], work_25);
      gsl_rstat_quantile_add(data[i], work_50);
      gsl_rstat_quantile_add(data[i], work_75);
    }

  /* exact values */
  gsl_sort(data, 1, N);
  exact_p25 = gsl_stats_quantile_from_sorted_data(data, 1, N, 0.25);
  exact_p50 = gsl_stats_quantile_from_sorted_data(data, 1, N, 0.5);
  exact_p75 = gsl_stats_quantile_from_sorted_data(data, 1, N, 0.75);

  /* estimated values */
  val_p25 = gsl_rstat_quantile_get(work_25);
  val_p50 = gsl_rstat_quantile_get(work_50);
  val_p75 = gsl_rstat_quantile_get(work_75);

  printf ("The dataset is %g, %g, %g, %g, %g, ...\n",
         data[0], data[1], data[2], data[3], data[4]);

  printf ("0.25 quartile: exact = %.5f, estimated = %.5f, error = %.6e\n",
          exact_p25, val_p25, (val_p25 - exact_p25) / exact_p25);
  printf ("0.50 quartile: exact = %.5f, estimated = %.5f, error = %.6e\n",
          exact_p50, val_p50, (val_p50 - exact_p50) / exact_p50);
  printf ("0.75 quartile: exact = %.5f, estimated = %.5f, error = %.6e\n",
          exact_p75, val_p75, (val_p75 - exact_p75) / exact_p75);

  gsl_rstat_quantile_free(work_25);
  gsl_rstat_quantile_free(work_50);
  gsl_rstat_quantile_free(work_75);
  gsl_rng_free(r);
  free(data);

  return 0;
}
