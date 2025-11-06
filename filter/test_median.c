/* filter/test_median.c
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
#include <gsl/gsl_filter.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_statistics.h>

/* compute filtered data by explicitely constructing window, sorting it and finding median */
int
slow_movmedian(const gsl_filter_end_t etype, const gsl_vector * x, gsl_vector * y, const int K)
{
  const size_t H = K / 2;
  const size_t n = x->size;
  double *window = malloc(K * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill((gsl_movstat_end_t) etype, x, i, H, H, window);
      double yi = gsl_stats_median(window, 1, wsize);

      gsl_vector_set(y, i, yi);
    }

  free(window);

  return GSL_SUCCESS;
}

static void
test_median_proc(const double tol, const size_t n, const size_t K,
                 const gsl_filter_end_t etype, gsl_rng *rng_p)
{
  gsl_filter_median_workspace *w = gsl_filter_median_alloc(K);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *z = gsl_vector_alloc(n);
  char buf[2048];

  /* test moving median with random input */
  random_vector(x, rng_p);

  /* y = median(x) with slow brute force algorithm */
  slow_movmedian(etype, x, y, K);

  /* z = median(x) */
  gsl_filter_median(etype, x, z, w);

  /* test y = z */
  sprintf(buf, "n=%zu K=%zu endtype=%u median random", n, K, etype);
  compare_vectors(tol, z, y, buf);

  /* z = median(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_filter_median(etype, z, z, w);

  sprintf(buf, "n=%zu K=%zu endtype=%u median random in-place", n, K, etype);
  compare_vectors(tol, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_filter_median_free(w);
}

static void
test_median(gsl_rng * rng_p)
{
  test_median_proc(GSL_DBL_EPSILON, 1000, 1, GSL_FILTER_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 3, GSL_FILTER_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 5, GSL_FILTER_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 100, 301, GSL_FILTER_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 17, GSL_FILTER_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 9, GSL_FILTER_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 901, GSL_FILTER_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 201, GSL_FILTER_END_PADZERO, rng_p);

  test_median_proc(GSL_DBL_EPSILON, 1000, 1, GSL_FILTER_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 3, GSL_FILTER_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 5, GSL_FILTER_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 100, 301, GSL_FILTER_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 17, GSL_FILTER_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 9, GSL_FILTER_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 901, GSL_FILTER_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 201, GSL_FILTER_END_PADVALUE, rng_p);

  test_median_proc(GSL_DBL_EPSILON, 1000, 1, GSL_FILTER_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 3, GSL_FILTER_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 5, GSL_FILTER_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 100, 301, GSL_FILTER_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 17, GSL_FILTER_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 9, GSL_FILTER_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 901, GSL_FILTER_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 201, GSL_FILTER_END_TRUNCATE, rng_p);
}
