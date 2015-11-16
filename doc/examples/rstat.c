#include <stdio.h>
#include <gsl/gsl_rstat.h>

int
main(void)
{
  double data[5] = {17.2, 18.1, 16.5, 18.3, 12.6};
  double mean, variance, largest, smallest;
  gsl_rstat_workspace *rstat_p = gsl_rstat_alloc();
  size_t i;

  /* add data to rstat accumulator */
  for (i = 0; i < 5; ++i)
    gsl_rstat_add(data[i], rstat_p);

  mean     = gsl_rstat_mean(rstat_p);
  variance = gsl_rstat_variance(rstat_p);
  largest  = gsl_rstat_max(rstat_p);
  smallest = gsl_rstat_min(rstat_p);

  printf ("The dataset is %g, %g, %g, %g, %g\n",
         data[0], data[1], data[2], data[3], data[4]);

  printf ("The sample mean is %g\n", mean);
  printf ("The estimated variance is %g\n", variance);
  printf ("The largest value is %g\n", largest);
  printf ("The smallest value is %g\n", smallest);

  gsl_rstat_free(rstat_p);

  return 0;
}
