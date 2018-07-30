/* movstat/test.c
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
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>

/* compare two vectors */
static void
compare_vectors(const double tol, const gsl_vector * v, const gsl_vector * expected,
                const char * desc)
{
  const size_t n = v->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double vi = gsl_vector_get(v, i);
      double ui = gsl_vector_get(expected, i);

      gsl_test_rel(vi, ui, tol, "%s i=%zu", desc, i);
    }
}

/* generate random vector with elements in [-1,1] */
static void
random_vector(gsl_vector * v, gsl_rng * r)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      double vi = 2.0 * gsl_rng_uniform(r) - 1.0; /* in [-1,1] */
      gsl_vector_set(v, i, vi);
    }
}

static int
test_noisy_sine(const double sigma, gsl_vector * x, gsl_rng * r)
{
  const size_t n = x->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double ti = (double) i / (n - 1.0);
      double xi = sin(2.0 * M_PI * ti) +
                  sin(2.0 * M_PI * 40.0 * ti);
      double ei = gsl_ran_gaussian(r, sigma);

      gsl_vector_set(x, i, xi + ei);
    }

  return GSL_SUCCESS;
}

#include "test_mad.c"
#include "test_mean.c"
#include "test_median.c"
#include "test_minmax.c"
#include "test_Qn.c"
#include "test_qqr.c"
#include "test_sum.c"
#include "test_Sn.c"
#include "test_variance.c"

int
main()
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  test_mean(r);
  test_median(r);
  test_minmax(r);
  test_mad(r);
  test_Qn(r);
  test_qqr(r);
  test_sum(r);
  test_Sn(r);
  test_variance(r);

  gsl_rng_free(r);

  exit (gsl_test_summary());
}
