/* filter/test_gaussian.c
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

/* compute Gaussian filter by explicitely constructing window and computing weighted sum */
int
slow_gaussian(const gsl_filter_end_t etype, const double alpha, const size_t order, const gsl_vector * x,
              gsl_vector * y, const size_t K)
{
  const size_t n = x->size;
  const size_t H = K / 2;
  double *window = malloc(K * sizeof(double));
  double *kernel = malloc(K * sizeof(double));
  gsl_vector_view k = gsl_vector_view_array(kernel, K);
  size_t i;

  gsl_filter_gaussian_kernel(alpha, order, 1, &k.vector);

  for (i = 0; i < n; ++i)
    {
      size_t wsize = gsl_movstat_fill((gsl_movstat_end_t) etype, x, i, H, H, window);
      double sum = 0.0;
      size_t j;

      for (j = 0; j < wsize; ++j)
        sum += window[j] * kernel[wsize - j - 1];

      gsl_vector_set(y, i, sum);
    }

  free(window);
  free(kernel);

  return GSL_SUCCESS;
}

static void
fdiff(const gsl_vector * x, gsl_vector * dx)
{
  const size_t N = x->size;
  size_t i;

  for (i = 1; i < N - 1; ++i)
    {
      double xm1 = gsl_vector_get(x, i - 1);
      double xp1 = gsl_vector_get(x, i + 1);
      gsl_vector_set(dx, i, 0.5 * (xp1 - xm1));
    }

  gsl_vector_set(dx, 0, gsl_vector_get(x, 1) - gsl_vector_get(x, 0));
  gsl_vector_set(dx, N - 1, gsl_vector_get(x, N - 1) - gsl_vector_get(x, N - 2));
}

static void
test_gaussian_kernel(const double alpha, const size_t K)
{
  const size_t max_order = 3;
  gsl_vector * kernel = gsl_vector_alloc(K);
  gsl_vector * deriv = gsl_vector_alloc(K);
  gsl_vector * deriv_fd = gsl_vector_alloc(K);
  char buf[2048];
  size_t order;

  gsl_filter_gaussian_kernel(alpha, 0, 0, kernel);

  for (order = 1; order <= max_order; ++order)
    {
      gsl_filter_gaussian_kernel(alpha, order, 0, deriv);
      fdiff(kernel, deriv_fd);

      sprintf(buf, "gaussian kernel order=%zu alpha=%g K=%zu", order, alpha, K);
      compare_vectors(1.0e-2, deriv_fd, deriv, buf);

      gsl_vector_memcpy(kernel, deriv);
    }

  gsl_vector_free(kernel);
  gsl_vector_free(deriv);
  gsl_vector_free(deriv_fd);
}

static void
test_gaussian_proc(const double tol, const double alpha, const size_t order, const size_t n, const size_t K,
                   const gsl_filter_end_t etype, gsl_rng * rng_p)
{
  gsl_filter_gaussian_workspace * w = gsl_filter_gaussian_alloc(K);
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * y = gsl_vector_alloc(n);
  gsl_vector * z = gsl_vector_alloc(n);
  char buf[2048];

  random_vector(x, rng_p);

  /* y = filter(x) with slow brute force method */
  slow_gaussian(etype, alpha, order, x, y, K);

  /* y = filter(x) with fast method */
  gsl_filter_gaussian(etype, alpha, order, x, z, w);

  /* test y = z */
  sprintf(buf, "n=%zu K=%zu endtype=%u alpha=%g order=%zu gaussian random", n, K, etype, alpha, order);
  compare_vectors(tol, z, y, buf);

  /* z = filter(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_filter_gaussian(etype, alpha, order, z, z, w);

  sprintf(buf, "n=%zu K=%zu endtype=%u alpha=%g order=%zu gaussian random in-place", n, K, etype, alpha, order);
  compare_vectors(tol, z, y, buf);

  gsl_filter_gaussian_free(w);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
}

#if 0
static void
test_gaussian_deriv(const double alpha, const size_t n, const size_t K)
{
  const double f_low = 1.0;
  const double f_high = 50.0;
  const double gamma = 2.0 * M_PI / (n - 1.0);
  const double dt = 1.0 / (n - 1.0);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *dx = gsl_vector_alloc(n);
  gsl_vector *y1 = gsl_vector_alloc(n);
  gsl_vector *y2 = gsl_vector_alloc(n);
  gsl_filter_gaussian_workspace *w = gsl_filter_gaussian_alloc(K);
  size_t i;

  /* make input signal composed of two sine waves at different frequencies */
  for (i = 0; i < n; ++i)
    {
      double xi = sin(gamma * f_low * i) + sin(gamma * f_high * i);
      double dxi = gamma * f_low * cos(gamma * f_low * i) +
                   gamma * f_high * cos(gamma * f_high * i);

      gsl_vector_set(x, i, xi);
      gsl_vector_set(dx, i, dxi);
    }

  /* compute y1 = G * dx(t)/dt */
  gsl_filter_gaussian(alpha, 0, dx, y1, w);

  /* compute y2 = dG/dt * x(t) */
  gsl_filter_gaussian(alpha, 1, x, y2, w);

  for (i = 0; i < n; ++i)
    {
      printf("%zu %.12e %.12e %.12e %.12e\n",
             i,
             gsl_vector_get(x, i),
             gsl_vector_get(dx, i),
             gsl_vector_get(y1, i),
             gsl_vector_get(y2, i));
    }

  gsl_vector_free(x);
  gsl_vector_free(dx);
  gsl_vector_free(y1);
  gsl_vector_free(y2);
  gsl_filter_gaussian_free(w);
}
#endif

static void
test_gaussian(gsl_rng * r)
{
  const double tol = 1.0e-10;
  size_t order;

  test_gaussian_kernel(3.0, 2001);

  for (order = 0; order <= 3; ++order)
    {
      test_gaussian_proc(tol, 2.5, order, 1000, 21, GSL_FILTER_END_PADZERO, r);
      test_gaussian_proc(tol, 3.0, order, 500, 11, GSL_FILTER_END_PADZERO, r);
      test_gaussian_proc(tol, 1.0, order, 50, 101, GSL_FILTER_END_PADZERO, r);
      test_gaussian_proc(tol, 2.0, order, 50, 11, GSL_FILTER_END_PADZERO, r);

      test_gaussian_proc(tol, 2.5, order, 1000, 21, GSL_FILTER_END_PADVALUE, r);
      test_gaussian_proc(tol, 3.0, order, 500, 11, GSL_FILTER_END_PADVALUE, r);
      test_gaussian_proc(tol, 1.0, order, 50, 101, GSL_FILTER_END_PADVALUE, r);
      test_gaussian_proc(tol, 2.0, order, 50, 11, GSL_FILTER_END_PADVALUE, r);

      test_gaussian_proc(tol, 2.5, order, 1000, 21, GSL_FILTER_END_TRUNCATE, r);
      test_gaussian_proc(tol, 3.0, order, 500, 11, GSL_FILTER_END_TRUNCATE, r);
      test_gaussian_proc(tol, 1.0, order, 50, 101, GSL_FILTER_END_TRUNCATE, r);
      test_gaussian_proc(tol, 2.0, order, 50, 11, GSL_FILTER_END_TRUNCATE, r);
    }
}
