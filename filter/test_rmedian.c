/* filter/test_rmedian.c
 * 
 * Copyright (C) 2018 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

/* slow/dumb rmedian which constructs actual window for each sample, sorts
 * it and finds median */
static int
slow_rmedian(const gsl_filter_end_t endtype, const gsl_vector * x, gsl_vector * y, const int K)
{
  const int n = (int) x->size;
  const int H = K / 2;
  double *window = malloc(K * sizeof(double));
  int i;

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill((gsl_movstat_end_t) endtype, x, i, H, H, window);
      double yi;
      int j;

      /* fill first half of window with previous filter output values */
      for (j = i - H; j < i; ++j)
        {
          if (j < 0)
            {
              if (endtype == GSL_FILTER_END_PADVALUE)
                window[j - i + H] = gsl_vector_get(x, 0);
              else if (endtype == GSL_FILTER_END_PADZERO)
                window[j - i + H] = 0.0;
            }
          else
            {
              window[j - i + H] = gsl_vector_get(y, j);
            }
        }

      yi = gsl_stats_median(window, 1, wsize);
      gsl_vector_set(y, i, yi);
    }

  free(window);

  return GSL_SUCCESS;
}

/* test square wave input (root signal) */
static void
test_rmedian_root(const gsl_filter_end_t etype, const size_t n, const size_t k)
{
  const double tol = 1.0e-12;
  gsl_filter_rmedian_workspace *w = gsl_filter_rmedian_alloc(k);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  char buf[2048];
  size_t i;

  /* test a root sequence (square input): x = [zero one zero] */
  gsl_vector_set_all(x, 0.0);

  for (i = n / 3; i <= n / 2; ++i)
    gsl_vector_set(x, i, 1.0);

  /* compute y = rmedian(x) and test y = x */
  gsl_filter_rmedian(etype, x, y, w);

  sprintf(buf, "n=%zu k=%zu RMF square wave root sequence", n, k);
  compare_vectors(tol, y, x, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_filter_rmedian_free(w);
}

/* test random input and in-place */
static void
test_rmedian_random(const gsl_filter_end_t etype, const size_t n, const int K, gsl_rng * r)
{
  const double tol = 1.0e-12;
  gsl_filter_rmedian_workspace *w = gsl_filter_rmedian_alloc(K);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *z = gsl_vector_alloc(n);
  char buf[2048];

  /* test filter with random input against slow algorithm */
  random_vector(x, r);

  /* y = rmedian(x) */
  gsl_filter_rmedian(etype, x, y, w);

  /* y = rmedian(x) with slow algorithm */
  slow_rmedian(etype, x, z, K);

  /* test y = z */
  sprintf(buf, "n=%zu K=%d RMF symmetric random slow test", n, K);
  compare_vectors(tol, y, z, buf);

  /* test in-place filter */

  /* z = rmedian(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_filter_rmedian(etype, z, z, w);

  sprintf(buf, "n=%zu K=%d RMF symmetric random in-place", n, K);
  compare_vectors(tol, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_filter_rmedian_free(w);
}

void
test_rmedian(gsl_rng * rng_p)
{
  /* test root sequences */

  test_rmedian_root(GSL_FILTER_END_PADZERO, 1000, 3);
  test_rmedian_root(GSL_FILTER_END_PADZERO, 2000, 101);

  test_rmedian_root(GSL_FILTER_END_PADVALUE, 1000, 3);
  test_rmedian_root(GSL_FILTER_END_PADVALUE, 2000, 101);

  /* test random input */

  test_rmedian_random(GSL_FILTER_END_PADZERO, 10, 1, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADZERO, 100, 3, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADZERO, 1000, 3, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADZERO, 100, 1001, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADZERO, 5, 7, rng_p);

  test_rmedian_random(GSL_FILTER_END_PADVALUE, 10, 1, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADVALUE, 100, 3, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADVALUE, 1000, 3, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADVALUE, 100, 1001, rng_p);
  test_rmedian_random(GSL_FILTER_END_PADVALUE, 5, 7, rng_p);
}
