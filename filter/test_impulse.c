/* filter/test_impulse.c
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
#include <gsl/gsl_filter.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static int
vector_sum(const gsl_vector_int * v)
{
  size_t i;
  int sum = 0;

  for (i = 0; i < v->size; ++i)
    {
      int vi = gsl_vector_int_get(v, i);
      sum += vi;
    }

  return sum;
}

static void
test_impulse_proc(const double tol, const size_t n, const size_t K, const double nsigma,
                  const gsl_filter_end_t etype, const gsl_filter_scale_t stype,
                  const double outlier_percentage, gsl_rng * r)
{
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * y = gsl_vector_alloc(n);
  gsl_vector * z = gsl_vector_alloc(n);
  gsl_vector * y_med = gsl_vector_alloc(n);
  gsl_vector * xmedian = gsl_vector_alloc(n);
  gsl_vector * xsigma = gsl_vector_alloc(n);
  size_t noutlier;
  gsl_vector_int * ioutlier = gsl_vector_int_alloc(n);
  gsl_vector_int * ioutlier_exact = gsl_vector_int_alloc(n);
  size_t i;
  gsl_filter_impulse_workspace *impulse_p = gsl_filter_impulse_alloc(K);
  gsl_filter_median_workspace *median_p = gsl_filter_median_alloc(K);
  size_t noutlier_exact = 0;
  char buf[1024];

  gsl_vector_int_set_zero(ioutlier_exact);

  for (i = 0; i < n; ++i)
    {
      double xi = gsl_ran_gaussian(r, 1.0);
      double vi = gsl_rng_uniform(r);

      if (vi <= outlier_percentage)
        {
          xi += 15.0 * GSL_SIGN(xi);
          ++noutlier_exact;
          gsl_vector_int_set(ioutlier_exact, i, 1);
        }

      gsl_vector_set(x, i, xi);
    }

  /* first test that median filter is equal to impulse filter with nsigma = 0 */

  gsl_filter_median(etype, x, y_med, median_p);
  gsl_filter_impulse(etype, stype, 0.0, x, y, xmedian, xsigma, &noutlier, ioutlier, impulse_p);

  sprintf(buf, "impulse nsigma=0 smf comparison, etype=%u stype=%u", etype, stype);
  compare_vectors(tol, y, y_med, buf);

  /* second test: filter y = impulse(x) with given nsigma */

  gsl_filter_impulse(etype, stype, nsigma, x, y, xmedian, xsigma, &noutlier, ioutlier, impulse_p);

  /* test correct number of outliers detected */
  gsl_test(noutlier != noutlier_exact, "impulse [n=%zu,K=%zu,nsigma=%g,outlier_percentage=%g] noutlier=%zu exact=%zu",
           n, K, nsigma, outlier_percentage, noutlier, noutlier_exact);

#if 0
  {
    for (i = 0; i < n; ++i)
      {
        printf("%.12e %.12e %d %.12e %.12e\n",
               gsl_vector_get(x, i),
               gsl_vector_get(y, i),
               gsl_vector_int_get(ioutlier, i),
               gsl_vector_get(xmedian, i) + nsigma * gsl_vector_get(xsigma, i),
               gsl_vector_get(xmedian, i) - nsigma * gsl_vector_get(xsigma, i));
      }
  }
#endif

  /* test outliers found in correct locations */
  for (i = 0; i < n; ++i)
    {
      int val = gsl_vector_int_get(ioutlier, i);
      int val_exact = gsl_vector_int_get(ioutlier_exact, i);

      gsl_test(val != val_exact, "test_impulse: outlier mismatch [i=%zu,K=%zu,nsigma=%g,outlier_percentage=%g] ioutlier=%d ioutlier_exact=%d",
               i, K, nsigma, outlier_percentage, val, val_exact);
    }

  /* test noutlier = sum(ioutlier) */
  {
    size_t iout_sum = vector_sum(ioutlier);
    gsl_test(noutlier != iout_sum, "impulse [K=%zu,nsigma=%g,outlier_percentage=%g] noutlier=%zu sum(ioutlier)=%zu",
             K, nsigma, outlier_percentage, noutlier, iout_sum);
  }

  /* third test: test in-place filter z = impulse(z) */

  gsl_vector_memcpy(z, x);
  gsl_filter_impulse(etype, stype, nsigma, z, z, xmedian, xsigma, &noutlier, ioutlier, impulse_p);

  sprintf(buf, "impulse in-place nsigma=%g,n=%zu,K=%zu,etype=%u stype=%u", nsigma, n, K, etype, stype);
  compare_vectors(GSL_DBL_EPSILON, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_vector_free(y_med);
  gsl_vector_free(xmedian);
  gsl_vector_free(xsigma);
  gsl_vector_int_free(ioutlier);
  gsl_vector_int_free(ioutlier_exact);
  gsl_filter_impulse_free(impulse_p);
  gsl_filter_median_free(median_p);
}

static void
test_impulse(gsl_rng * r)
{
  const double tol = 1.0e-10;

  test_impulse_proc(tol, 1000, 21, 6.0, GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_QN, 0.05, r);
  test_impulse_proc(tol, 1000, 21, 6.0, GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_SN, 0.05, r);
  test_impulse_proc(tol, 1000, 21, 6.0, GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_MAD, 0.05, r);
  test_impulse_proc(tol, 1000, 21, 6.0, GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_IQR, 0.05, r);
}
