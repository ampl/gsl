/* movstat/test_variance.c
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

/* compute moving variance by explicitely constructing window and computing variance */
int
slow_movvar(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y,
            const int H, const int J)
{
  const size_t n = x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill(etype, x, i, H, J, window);
      double variance = (wsize > 1) ? gsl_stats_variance(window, 1, wsize) : 0.0;

      gsl_vector_set(y, i, variance);
    }

  free(window);

  return GSL_SUCCESS;
}

/* compute moving variance by explicitely constructing window and computing variance */
int
slow_movsd(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y,
           const int H, const int J)
{
  const size_t n = x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill(etype, x, i, H, J, window);
      double sd = (wsize > 1) ? gsl_stats_sd(window, 1, wsize) : 0.0;

      gsl_vector_set(y, i, sd);
    }

  free(window);

  return GSL_SUCCESS;
}

static double
func_var(const size_t n, double x[], void * params)
{
  (void) params;
  if (n > 1)
    return gsl_stats_variance(x, 1, n);
  else
    return 0.0;
}

static double
func_sd(const size_t n, double x[], void * params)
{
  (void) params;
  if (n > 1)
    return gsl_stats_sd(x, 1, n);
  else
    return 0.0;
}

static void
test_variance_proc(const double tol, const size_t n, const size_t H, const size_t J,
                   const gsl_movstat_end_t etype, gsl_rng * rng_p)
{
  gsl_movstat_workspace * w = gsl_movstat_alloc2(H, J);
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * y = gsl_vector_alloc(n);
  gsl_vector * z = gsl_vector_alloc(n);
  gsl_movstat_function F1, F2;
  char buf[2048];

  F1.function = func_var;
  F1.params = NULL;

  F2.function = func_sd;
  F2.params = NULL;

  random_vector(x, rng_p);

  /* test variance */

  /* y = variance(x) with slow brute force method */
  slow_movvar(etype, x, y, H, J);

  /* y = variance(x) with fast method */
  gsl_movstat_variance(etype, x, z, w);

  /* test y = z */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u variance random", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = variance(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_variance(etype, z, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u variance random in-place", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = variance(x) with user-defined function */
  gsl_movstat_apply(etype, &F1, x, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u variance user", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* test standard deviation */

  /* y = stddev(x) with slow brute force method */
  slow_movsd(etype, x, y, H, J);

  /* y = stddev(x) with fast method */
  gsl_movstat_sd(etype, x, z, w);

  /* test y = z */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u stddev random", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = stddev(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_sd(etype, z, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u stddev random in-place", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = stddev(x) with user-defined function */
  gsl_movstat_apply(etype, &F2, x, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u stddev user", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  gsl_movstat_free(w);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
}

static void
test_variance(gsl_rng * rng_p)
{
  const double eps = 1.0e-10;

  test_variance_proc(eps, 1000, 0, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_variance_proc(eps, 1000, 3, 3, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_variance_proc(eps, 1000, 0, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_variance_proc(eps, 1000, 5, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_variance_proc(eps, 2000, 10, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_variance_proc(eps, 2000, 5, 10, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_variance_proc(eps, 20, 50, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_variance_proc(eps, 20, 10, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_variance_proc(eps, 20, 50, 10, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_variance_proc(eps, 1000, 0, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_variance_proc(eps, 1000, 3, 3, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_variance_proc(eps, 1000, 1, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_variance_proc(eps, 1000, 5, 1, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_variance_proc(eps, 2000, 10, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_variance_proc(eps, 2000, 5, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_variance_proc(eps, 20, 50, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_variance_proc(eps, 20, 10, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_variance_proc(eps, 20, 50, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_variance_proc(eps, 1000, 0, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_variance_proc(eps, 1000, 3, 3, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_variance_proc(eps, 1000, 0, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_variance_proc(eps, 1000, 5, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_variance_proc(eps, 2000, 10, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_variance_proc(eps, 2000, 5, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_variance_proc(eps, 20, 50, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_variance_proc(eps, 20, 10, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_variance_proc(eps, 20, 50, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
}
