/* movstat/test_median.c
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
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_statistics.h>

/* compute filtered data by explicitely constructing window, sorting it and finding median */
int
slow_movmedian(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y,
               const int H, const int J)
{
  const size_t n = x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill(etype, x, i, H, J, window);
      double yi = gsl_stats_median(window, 1, wsize);

      gsl_vector_set(y, i, yi);
    }

  free(window);

  return GSL_SUCCESS;
}

/* test root sequence */
static void
test_median_root(const double tol, const size_t n, const size_t K, const gsl_movstat_end_t etype)
{
  gsl_movstat_workspace *w = gsl_movstat_alloc(K);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  char buf[2048];
  size_t i;

  /* test a root sequence (square input): x = [zero one zero] */
  gsl_vector_set_all(x, 0.0);

  for (i = n / 3; i <= n / 2; ++i)
    gsl_vector_set(x, i, 1.0);

  /* compute y = median(x) and test y = x */
  gsl_movstat_median(etype, x, y, w);

  sprintf(buf, "n=%zu K=%zu endtype=%u SMF root sequence", n, K, etype);
  compare_vectors(tol, y, x, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_movstat_free(w);
}

static double
func_median(const size_t n, double x[], void * params)
{
  (void) params;
  return gsl_stats_median(x, 1, n);
}

static void
test_median_proc(const double tol, const size_t n, const size_t H, const size_t J,
                 const gsl_movstat_end_t etype, gsl_rng *rng_p)
{
  gsl_movstat_workspace *w;
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *z = gsl_vector_alloc(n);
  gsl_movstat_function F;
  char buf[2048];

  F.function = func_median;
  F.params = NULL;

  if (H == J)
    w = gsl_movstat_alloc(2*H + 1);
  else
    w = gsl_movstat_alloc2(H, J);

  /* test moving median with random input */
  random_vector(x, rng_p);

  /* y = median(x) with slow brute force algorithm */
  slow_movmedian(etype, x, y, H, J);

  /* z = median(x) */
  gsl_movstat_median(etype, x, z, w);

  /* test y = z */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u median random", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = median(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_median(etype, z, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u median random in-place", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = median(x) with user-defined function */
  gsl_movstat_apply(etype, &F, x, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u median user", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_movstat_free(w);
}

static void
test_median(gsl_rng * rng_p)
{
  test_median_root(GSL_DBL_EPSILON, 1000, 3, GSL_MOVSTAT_END_PADZERO);
  test_median_root(GSL_DBL_EPSILON, 200, 15, GSL_MOVSTAT_END_PADVALUE);
  test_median_root(GSL_DBL_EPSILON, 100, 5, GSL_MOVSTAT_END_TRUNCATE);

  test_median_proc(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 1, 1, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 100, 150, 150, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 8, 8, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 0, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 5, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 15, 10, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 10, 15, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 100, 150, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 150, 100, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 100, 100, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_median_proc(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 1, 1, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 100, 150, 150, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 8, 8, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 0, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 5, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 15, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 10, 15, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 100, 150, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 150, 100, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 100, 100, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_median_proc(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 1000, 1, 1, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 100, 150, 150, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 8, 8, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 0, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 5, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 15, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 5000, 10, 15, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 100, 150, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 150, 100, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_proc(GSL_DBL_EPSILON, 50, 100, 100, GSL_MOVSTAT_END_TRUNCATE, rng_p);
}
