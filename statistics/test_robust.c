/* statistics/test_robust.c
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

#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_ieee_utils.h>

int test_robust (void);

/* random vector in [-1,1] */
static int
random_array(const size_t n, double * x, gsl_rng * r)
{
  size_t i;

  for (i = 0; i < n; ++i)
    x[i] = 2.0 * gsl_rng_uniform(r) - 1.0;

  return 0;
}

/* calculate MAD statistic for input vector using slow/naive algorithm */
static double
slow_MAD(const size_t n, const double x[])
{
  double *work = malloc(n * sizeof(double));
  double median, mad;
  size_t i;

  for (i = 0; i < n; ++i)
    work[i] = x[i];

  gsl_sort(work, 1, n);
  median = gsl_stats_median_from_sorted_data(work, 1, n);

  for (i = 0; i < n; ++i)
    work[i] = fabs(x[i] - median);

  gsl_sort(work, 1, n);
  mad = gsl_stats_median_from_sorted_data(work, 1, n);

  free(work);

  return mad;
}

/* calculate S_n statistic for input vector using slow/naive algorithm */
static double
slow_Sn0(const size_t n, const double x[])
{
  double *work1 = malloc(n * sizeof(double));
  double *work2 = malloc(n * sizeof(double));
  double Sn;
  size_t i, j;

  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
        work1[j] = fabs(x[i] - x[j]);

      /* find himed_j | x_i - x_j | */
      gsl_sort(work1, 1, n);
      work2[i] = work1[n / 2];
    }

  /* find lomed_i { himed_j | x_i - x_j | } */
  gsl_sort(work2, 1, n);
  Sn = work2[(n + 1) / 2 - 1];

  free(work1);
  free(work2);

  return Sn;
}

/* calculate Q_n statistic for input vector using slow/naive algorithm */
static double
slow_Qn0(const size_t n, const double x[])
{
  const size_t wsize = n * (n - 1) / 2;
  const size_t n_2 = n / 2;
  const size_t k = ((n_2 + 1) * n_2) / 2;
  double *work;
  double Qn;
  size_t idx = 0;
  size_t i, j;

  if (n < 2)
    return (0.0);

  work = malloc(wsize * sizeof(double));

  for (i = 0; i < n; ++i)
    {
      for (j = i + 1; j < n; ++j)
        work[idx++] = fabs(x[i] - x[j]);
    }

  gsl_sort(work, 1, idx);
  Qn = work[k - 1];

  free(work);

  return Qn;
}

static int
test_median(const double tol, const size_t n, gsl_rng * r)
{
  double * x = malloc(n * sizeof(double));
  double median1, median2;

  random_array(n, x, r);

  median1 = gsl_stats_median(x, 1, n);

  gsl_sort(x, 1, n);
  median2 = gsl_stats_median_from_sorted_data(x, 1, n);

  gsl_test_rel(median1, median2, tol, "test_median n=%zu", n);

  free(x);

  return 0;
}

static int
test_mad(const double tol, const size_t n, gsl_rng * r)
{
  double * x = malloc(n * sizeof(double));
  double * work = malloc(n * sizeof(double));
  double mad1, mad2;

  random_array(n, x, r);

  mad1 = slow_MAD(n, x);

  gsl_sort(x, 1, n);
  mad2 = gsl_stats_mad0(x, 1, n, work);

  gsl_test_rel(mad1, mad2, tol, "test_mad n=%zu", n);

  free(x);
  free(work);

  return 0;
}

static int
test_Sn(const double tol, const size_t n, gsl_rng * r)
{
  double * x = malloc(n * sizeof(double));
  double * work = malloc(n * sizeof(double));
  double Sn1, Sn2;

  random_array(n, x, r);

  /* compute S_n with slow/naive algorithm */
  Sn1 = slow_Sn0(n, x);

  /* compute S_n with efficient algorithm */
  gsl_sort(x, 1, n);
  Sn2 = gsl_stats_Sn0_from_sorted_data(x, 1, n, work);

  gsl_test_rel(Sn2, Sn1, tol, "test_Sn n=%zu", n);

  free(x);
  free(work);

  return 0;
}

static int
test_Qn(const double tol, const size_t n, gsl_rng * r)
{
  double * x = malloc(n * sizeof(double));
  double * work = malloc(3 * n * sizeof(double));
  int * work_int = malloc(5 * n * sizeof(int));
  double Qn1, Qn2;

  random_array(n, x, r);

  /* compute Q_n with slow/naive algorithm */
  Qn1 = slow_Qn0(n, x);

  /* compute Q_n with efficient algorithm */
  gsl_sort(x, 1, n);
  Qn2 = gsl_stats_Qn0_from_sorted_data(x, 1, n, work, work_int);

  gsl_test_rel(Qn2, Qn1, tol, "test_Qn n=%zu", n);

  free(x);
  free(work);
  free(work_int);

  return 0;
}

int
test_robust (void)
{
  const double tol = 1.0e-12;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_default);

  test_median(GSL_DBL_EPSILON, 1, r);
  test_median(GSL_DBL_EPSILON, 2, r);
  test_median(GSL_DBL_EPSILON, 3, r);
  test_median(GSL_DBL_EPSILON, 100, r);
  test_median(GSL_DBL_EPSILON, 101, r);
  test_median(GSL_DBL_EPSILON, 500, r);
  test_median(GSL_DBL_EPSILON, 501, r);

  test_mad(GSL_DBL_EPSILON, 1, r);
  test_mad(GSL_DBL_EPSILON, 2, r);
  test_mad(GSL_DBL_EPSILON, 3, r);
  test_mad(GSL_DBL_EPSILON, 100, r);
  test_mad(GSL_DBL_EPSILON, 101, r);
  test_mad(GSL_DBL_EPSILON, 500, r);
  test_mad(GSL_DBL_EPSILON, 501, r);

  test_Sn(tol, 1, r);
  test_Sn(tol, 2, r);
  test_Sn(tol, 3, r);
  test_Sn(tol, 100, r);
  test_Sn(tol, 101, r);
  test_Sn(tol, 500, r);
  test_Sn(tol, 501, r);

  test_Qn(tol, 1, r);
  test_Qn(tol, 2, r);
  test_Qn(tol, 3, r);
  test_Qn(tol, 4, r);
  test_Qn(tol, 5, r);
  test_Qn(tol, 100, r);
  test_Qn(tol, 101, r);
  test_Qn(tol, 500, r);
  test_Qn(tol, 501, r);

  gsl_rng_free(r);

  return 0;
}
