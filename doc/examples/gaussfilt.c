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
  const size_t N = 500;                        /* length of time series */
  const size_t K = 51;                         /* window size */
  const double alpha[3] = { 0.5, 3.0, 10.0 };  /* alpha values */
  gsl_vector *x = gsl_vector_alloc(N);         /* input vector */
  gsl_vector *y1 = gsl_vector_alloc(N);        /* filtered output vector for alpha1 */
  gsl_vector *y2 = gsl_vector_alloc(N);        /* filtered output vector for alpha2 */
  gsl_vector *y3 = gsl_vector_alloc(N);        /* filtered output vector for alpha3 */
  gsl_vector *k1 = gsl_vector_alloc(K);        /* Gaussian kernel for alpha1 */
  gsl_vector *k2 = gsl_vector_alloc(K);        /* Gaussian kernel for alpha2 */
  gsl_vector *k3 = gsl_vector_alloc(K);        /* Gaussian kernel for alpha3 */
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_filter_gaussian_workspace *gauss_p = gsl_filter_gaussian_alloc(K);
  size_t i;
  double sum = 0.0;

  /* generate input signal */
  for (i = 0; i < N; ++i)
    {
      double ui = gsl_ran_gaussian(r, 1.0);
      sum += ui;
      gsl_vector_set(x, i, sum);
    }

  /* compute kernels without normalization */
  gsl_filter_gaussian_kernel(alpha[0], 0, 0, k1);
  gsl_filter_gaussian_kernel(alpha[1], 0, 0, k2);
  gsl_filter_gaussian_kernel(alpha[2], 0, 0, k3);

  /* apply filters */
  gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha[0], 0, x, y1, gauss_p);
  gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha[1], 0, x, y2, gauss_p);
  gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha[2], 0, x, y3, gauss_p);

  /* print kernels */
  for (i = 0; i < K; ++i)
    {
      double k1i = gsl_vector_get(k1, i);
      double k2i = gsl_vector_get(k2, i);
      double k3i = gsl_vector_get(k3, i);

      printf("%e %e %e\n", k1i, k2i, k3i);
    }

  printf("\n\n");

  /* print filter results */
  for (i = 0; i < N; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double y1i = gsl_vector_get(y1, i);
      double y2i = gsl_vector_get(y2, i);
      double y3i = gsl_vector_get(y3, i);

      printf("%.12e %.12e %.12e %.12e\n", xi, y1i, y2i, y3i);
    }

  gsl_vector_free(x);
  gsl_vector_free(y1);
  gsl_vector_free(y2);
  gsl_vector_free(y3);
  gsl_vector_free(k1);
  gsl_vector_free(k2);
  gsl_vector_free(k3);
  gsl_rng_free(r);
  gsl_filter_gaussian_free(gauss_p);

  return 0;
}
