/* movstat/test_mad.c
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
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_statistics.h>

/* compute filtered data by explicitely constructing window and finding MAD */
int
slow_movmad(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y,
            const int H, const int J)
{
  const size_t n = x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill(etype, x, i, H, J, window);
      double median = gsl_stats_median(window, 1, wsize);
      double mad;
      size_t j;

      for (j = 0; j < wsize; ++j)
        window[j] = fabs(window[j] - median);

      mad = gsl_stats_median(window, 1, wsize);

      gsl_vector_set(y, i, mad);
    }

  free(window);

  return GSL_SUCCESS;
}

static double
func_mad(const size_t n, double x[], void * params)
{
  double *work = malloc(n * sizeof(double));
  double mad;

  (void) params;

  mad = gsl_stats_mad0(x, 1, n, work);

  free(work);

  return mad;
}

static void
test_mad_proc(const double tol, const size_t n, const size_t H, const size_t J,
              const gsl_movstat_end_t etype, gsl_rng *rng_p)
{
  const double scale = 1.482602218505602;
  gsl_movstat_workspace *w;
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *z = gsl_vector_alloc(n);
  gsl_vector *med1 = gsl_vector_alloc(n);
  gsl_vector *med2 = gsl_vector_alloc(n);
  gsl_movstat_function F;
  char buf[2048];

  F.function = func_mad;
  F.params = NULL;

  if (H == J)
    w = gsl_movstat_alloc(2*H + 1);
  else
    w = gsl_movstat_alloc2(H, J);

  /* test moving MAD with random input */
  random_vector(x, rng_p);

  /* y = MAD(x) with slow brute force algorithm */
  slow_movmad(etype, x, y, H, J);

  /* z = MAD(x) */
  gsl_movstat_mad0(etype, x, med1, z, w);

  /* med2 = median(x) */
  gsl_movstat_median(etype, x, med2, w);

  /* test y = z */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u MAD0 random", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* test med1 = med2 */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u MAD0 random median test", n, H, J, etype);
  compare_vectors(tol, med1, med2, buf);

  /* z = MAD(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_mad0(etype, z, med1, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u MAD0 random in-place", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = MAD(x) with user-defined function */
  gsl_movstat_apply(etype, &F, x, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u MAD0 user", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* test scaled MAD function */

  /* z = MAD(x) */
  gsl_movstat_mad(etype, x, med1, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u MAD random", n, H, J, etype);
  gsl_vector_scale(y, scale);
  compare_vectors(tol, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_vector_free(med1);
  gsl_vector_free(med2);
  gsl_movstat_free(w);
}

static void
test_mad(gsl_rng * rng_p)
{
  test_mad_proc(GSL_DBL_EPSILON, 100, 0, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 1000, 1, 1, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 8, 8, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 0, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 5, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 15, 10, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 10, 15, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 100, 150, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 150, 100, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 100, 100, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_mad_proc(GSL_DBL_EPSILON, 100, 0, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 1000, 1, 1, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 8, 8, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 0, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 5, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 15, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 10, 15, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 100, 150, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 150, 100, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 100, 100, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_mad_proc(GSL_DBL_EPSILON, 100, 0, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 1000, 1, 1, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 8, 8, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 0, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 5, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 15, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 5000, 10, 15, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 100, 150, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 150, 100, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_mad_proc(GSL_DBL_EPSILON, 50, 100, 100, GSL_MOVSTAT_END_TRUNCATE, rng_p);
}
