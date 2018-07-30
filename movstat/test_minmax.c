/* movstat/test_minmax.c
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

/* compute filtered data by explicitely constructing window and finding min/max */
int
slow_minmax(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y_min, gsl_vector * y_max,
            const int H, const int J)
{
  const size_t n = x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill(etype, x, i, H, J, window);
      gsl_vector_view v = gsl_vector_view_array(window, wsize);
      double min, max;

      gsl_vector_minmax(&v.vector, &min, &max);
      gsl_vector_set(y_min, i, min);
      gsl_vector_set(y_max, i, max);
    }

  free(window);

  return GSL_SUCCESS;
}

static double
func_min(const size_t n, double x[], void * params)
{
  gsl_vector_view v = gsl_vector_view_array(x, n);
  (void) params;
  return gsl_vector_min(&v.vector);
}

static double
func_max(const size_t n, double x[], void * params)
{
  gsl_vector_view v = gsl_vector_view_array(x, n);
  (void) params;
  return gsl_vector_max(&v.vector);
}

static void
test_minmax_x(const double tol, const gsl_vector * x, const int H, const int J,
              const gsl_movstat_end_t endtype, const char * desc)
{
  const size_t n = x->size;
  gsl_vector * u_min = gsl_vector_alloc(n);
  gsl_vector * y_min = gsl_vector_alloc(n);
  gsl_vector * z_min = gsl_vector_alloc(n);
  gsl_vector * u_max = gsl_vector_alloc(n);
  gsl_vector * y_max = gsl_vector_alloc(n);
  gsl_vector * z_max = gsl_vector_alloc(n);
  gsl_movstat_workspace * w = gsl_movstat_alloc2(H, J);
  gsl_movstat_function F1, F2;
  char buf[2048];

  F1.function = func_min;
  F1.params = NULL;

  F2.function = func_max;
  F2.params = NULL;

  /* compute moving min/max */
  gsl_movstat_min(endtype, x, u_min, w);
  gsl_movstat_max(endtype, x, u_max, w);
  gsl_movstat_minmax(endtype, x, y_min, y_max, w);

  /* compute moving min/max with slow brute force method */
  slow_minmax(endtype, x, z_min, z_max, H, J);

  sprintf(buf, "test_minmax: %s min endtype=%d n=%zu H=%d J=%d", desc, endtype, n, H, J);
  compare_vectors(tol, u_min, z_min, buf);

  sprintf(buf, "test_minmax: %s max endtype=%d n=%zu H=%d J=%d", desc, endtype, n, H, J);
  compare_vectors(tol, u_max, z_max, buf);

  sprintf(buf, "test_minmax: %s minmax(minimum) endtype=%d n=%zu H=%d J=%d", desc, endtype, n, H, J);
  compare_vectors(tol, y_min, z_min, buf);

  sprintf(buf, "test_minmax: %s minmax(maximum) endtype=%d n=%zu H=%d J=%d", desc, endtype, n, H, J);
  compare_vectors(tol, y_max, z_max, buf);

  /* in-place tests */
  
  gsl_vector_memcpy(u_min, x);
  gsl_vector_memcpy(u_max, x);

  gsl_movstat_min(endtype, u_min, u_min, w);
  gsl_movstat_max(endtype, u_max, u_max, w);

  sprintf(buf, "test_minmax: %s in-place min endtype=%d n=%zu H=%d J=%d", desc, endtype, n, H, J);
  compare_vectors(tol, u_min, z_min, buf);

  sprintf(buf, "test_minmax: %s in-place max endtype=%d n=%zu H=%d J=%d", desc, endtype, n, H, J);
  compare_vectors(tol, u_max, z_max, buf);

  /* user-defined function tests */

  gsl_movstat_apply(endtype, &F1, x, z_min, w);
  sprintf(buf, "n=%zu H=%d J=%d endtype=%u min user", n, H, J, endtype);
  compare_vectors(tol, z_min, y_min, buf);

  gsl_movstat_apply(endtype, &F2, x, z_max, w);
  sprintf(buf, "n=%zu H=%d J=%d endtype=%u max user", n, H, J, endtype);
  compare_vectors(tol, z_max, y_max, buf);

  gsl_vector_free(u_min);
  gsl_vector_free(y_min);
  gsl_vector_free(z_min);
  gsl_vector_free(u_max);
  gsl_vector_free(y_max);
  gsl_vector_free(z_max);
  gsl_movstat_free(w);
}

/* test alternating sequence [a,b,a,b,...] input */
static void
test_minmax_alt(const double tol, const size_t n, const int H, const int J,
                const gsl_movstat_end_t endtype)
{
  const double a = 5.0;
  const double b = -5.0;
  gsl_vector * x = gsl_vector_alloc(n);
  size_t i;

  for (i = 0; i < n; ++i)
    {
      if (i % 2 == 0)
        gsl_vector_set(x, i, a);
      else
        gsl_vector_set(x, i, b);
    }

  test_minmax_x(tol, x, H, J, endtype, "alternating");

  gsl_vector_free(x);
}

/* test noisy sine wave input */
static void
test_minmax_sine(const double tol, const size_t n, const int H, const int J,
                 const gsl_movstat_end_t endtype, gsl_rng * rng_p)
{
  gsl_vector * x = gsl_vector_alloc(n);

  /* construct noisy sine signal */
  test_noisy_sine(0.5, x, rng_p);

  test_minmax_x(tol, x, H, J, endtype, "noisy_sine");

  gsl_vector_free(x);
}

/* test random input */
static void
test_minmax_random(const double tol, const size_t n, const int H, const int J,
                   const gsl_movstat_end_t endtype, gsl_rng * rng_p)
{
  gsl_vector * x = gsl_vector_alloc(n);

  /* construct random input signal */
  random_vector(x, rng_p);

  test_minmax_x(tol, x, H, J, endtype, "random");

  gsl_vector_free(x);
}

static void
test_minmax(gsl_rng * rng_p)
{
  /* alternating input */

  test_minmax_alt(GSL_DBL_EPSILON, 1000, 7, 7, GSL_MOVSTAT_END_PADZERO);
  test_minmax_alt(GSL_DBL_EPSILON, 1000, 5, 2, GSL_MOVSTAT_END_PADZERO);
  test_minmax_alt(GSL_DBL_EPSILON, 500, 1, 3, GSL_MOVSTAT_END_PADZERO);
  test_minmax_alt(GSL_DBL_EPSILON, 20, 50, 10, GSL_MOVSTAT_END_PADZERO);
  test_minmax_alt(GSL_DBL_EPSILON, 20, 10, 50, GSL_MOVSTAT_END_PADZERO);

  /* noisy sine wave input */

  test_minmax_sine(GSL_DBL_EPSILON, 1000, 5, 7, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 2000, 0, 2, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 500, 3, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 10, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 50, 10, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_minmax_sine(GSL_DBL_EPSILON, 500, 5, 7, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 10, 20, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 3, 3, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 10, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 50, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_minmax_sine(GSL_DBL_EPSILON, 500, 5, 7, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 10, 20, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 3, 3, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 30, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 5, 30, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 10, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 20, 50, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);

  /* random input */

  test_minmax_random(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 1000, 5, 7, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 2000, 0, 2, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 3, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 10, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 5, 10, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 10, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 50, 10, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_minmax_random(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 1000, 5, 7, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 2000, 0, 2, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 3, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 10, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 5, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 10, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 50, 10, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_minmax_random(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 1000, 5, 7, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 2000, 0, 2, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 3, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 10, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 500, 5, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 10, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_random(GSL_DBL_EPSILON, 20, 50, 10, GSL_MOVSTAT_END_TRUNCATE, rng_p);
}
