/* movstat/test_Sn.c
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
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_test.h>

/* calculate S_n statistic for input vector using slow/naive algorithm */
static int
slow_movSn(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y,
           const int H, const int J)
{
  const size_t n = x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  double *work = malloc(K * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill(etype, x, i, H, J, window);
      double Sn;

      gsl_sort(window, 1, wsize);
      Sn = gsl_stats_Sn_from_sorted_data(window, 1, wsize, work);
      gsl_vector_set(y, i, Sn);
    }

  free(window);
  free(work);

  return GSL_SUCCESS;
}

static double
func_Sn(const size_t n, double x[], void * params)
{
  double *work = malloc(n * sizeof(double));
  double Sn;

  (void) params;

  gsl_sort(x, 1, n);
  Sn = gsl_stats_Sn_from_sorted_data(x, 1, n, work);

  free(work);
  
  return Sn;
}

static void
test_Sn_proc(const double tol, const size_t n, const size_t H, const size_t J,
             const gsl_movstat_end_t etype, gsl_rng *rng_p)
{
  gsl_movstat_workspace *w;
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *z = gsl_vector_alloc(n);
  gsl_movstat_function F;
  char buf[2048];

  F.function = func_Sn;
  F.params = NULL;

  if (H == J)
    w = gsl_movstat_alloc(2*H + 1);
  else
    w = gsl_movstat_alloc2(H, J);

  /* test moving median with random input */
  random_vector(x, rng_p);

  /* y = S_n(x) with slow brute force algorithm */
  slow_movSn(etype, x, y, H, J);

  /* z = S_n(x) */
  gsl_movstat_Sn(etype, x, z, w);

  /* test y = z */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u Sn random", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = S_n(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_Sn(etype, z, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u Sn random in-place", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = S_n(x) with user-defined function */
  gsl_movstat_apply(etype, &F, x, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u Sn user", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_movstat_free(w);
}

static void
test_Sn(gsl_rng * rng_p)
{
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 5, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 5, 2, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 2, 5, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 2000, 50, 0, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 2000, 0, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 1, 50, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 50, 1, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_Sn_proc(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 5, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 5, 2, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 2, 5, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 2000, 50, 0, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 2000, 0, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 1, 50, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 50, 1, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_Sn_proc(GSL_DBL_EPSILON, 1000, 0, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 5, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 5, 2, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 1000, 2, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 2000, 50, 0, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 2000, 0, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 50, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 1, 50, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_Sn_proc(GSL_DBL_EPSILON, 20, 50, 1, GSL_MOVSTAT_END_TRUNCATE, rng_p);
}
