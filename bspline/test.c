/* bspline/test.c
 *
 * Copyright (C) 2006, 2007, 2009 Brian Gough
 * Copyright (C) 2008, 2011 Rhys Ulerich
 * Copyright (C) 2018, 2019, 2020, 2021, 2022 Patrick Alken
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>

/* generate random vector with elements in [a,b] */
static void
random_vector(const double a, const double b, gsl_vector * v, gsl_rng * r)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      double vi = (b - a) * gsl_rng_uniform(r) + a;
      gsl_vector_set(v, i, vi);
    }
}

static void
test_bspline(gsl_bspline_workspace * w)
{
  gsl_vector *B;
  gsl_matrix *dB;
  size_t i, j;
  size_t n = 100;
  size_t ncoeffs = gsl_bspline_ncontrol(w);
  size_t order = gsl_bspline_order(w);
  size_t nbreak = gsl_bspline_nbreak(w);
  double a = gsl_bspline_breakpoint(0, w);
  double b = gsl_bspline_breakpoint(nbreak - 1, w);

  B  = gsl_vector_alloc(ncoeffs);
  dB = gsl_matrix_alloc(ncoeffs, 1);

  /* ensure B-splines form a partition of unity */
  for (i = 0; i < n; i++)
    {
      double xi = a + (b - a) * (i / (n - 1.0));
      double sum = 0;
      gsl_bspline_eval_basis(xi, B, w);

      for (j = 0; j < ncoeffs; j++)
        {
          double Bj = gsl_vector_get(B, j);
          int s = (Bj < 0 || Bj > 1);
          gsl_test(s, "basis-spline coefficient %u is in range [0,1] for x=%g", j, xi);
          sum += Bj;
        }

      gsl_test_rel(sum, 1.0, order * GSL_DBL_EPSILON,
                   "basis-spline order %u is normalized for x=%g", order, xi);
    }

  /* ensure B-splines 0th derivatives agree with regular evaluation */
  for (i = 0; i < n; i++)
    {
      double xi = a + (b - a) * (i / (n - 1.0));
      gsl_bspline_eval_basis(xi, B, w);
      gsl_bspline_eval_deriv_basis(xi, 0, dB, w);

      for (j = 0; j < ncoeffs; j++)
        {
          gsl_test_abs(gsl_matrix_get(dB, j, 0), gsl_vector_get(B, j),
                       GSL_DBL_EPSILON,
                       "b-spline order %d basis #%d evaluation and 0th derivative consistent for x=%g",
                       order, j, xi);
        }

    }

  gsl_vector_free(B);
  gsl_matrix_free(dB);
}

#include "test_gram.c"
#include "test_greville.c"
#include "test_interp.c"
#include "test_integ.c"
#include "test_ls.c"

int
main(void)
{
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);

  gsl_ieee_env_setup();

  {
    size_t order, breakpoints;

    for (order = 1; order < 10; order++)
      {
        for (breakpoints = 2; breakpoints < 100; breakpoints++)
          {
            double a = -1.23 * order, b = 45.6 * order;
            gsl_bspline_workspace *w = gsl_bspline_alloc(order, breakpoints);
            gsl_bspline_init_uniform(a, b, w);
            test_bspline(w);
            gsl_bspline_free(w);
          }
      }

    for (order = 1; order < 10; order++)
      {
        for (breakpoints = 2; breakpoints < 100; breakpoints++)
          {
            double a = -1.23 * order, b = 45.6 * order;
            gsl_bspline_workspace *w = gsl_bspline_alloc(order, breakpoints);
            gsl_vector *k = gsl_vector_alloc(breakpoints);
            size_t i;

            for (i = 0; i < breakpoints; i++)
              {
                double f, x;
                f = sqrt(i / (breakpoints - 1.0));
                x = (1 - f) * a + f * b;
                gsl_vector_set(k, i, x);
              }

            gsl_bspline_init_augment(k, w);
            test_bspline(w);
            gsl_vector_free(k);
            gsl_bspline_free(w);
          }
      }
  }

  /* Spot check known 0th, 1st, 2nd derivative
     evaluations for a particular k = 2 case.  */
  {
    size_t i, j; /* looping */

    const double xloc[4]     =  { 0.0,  1.0,  6.0,  7.0};
    const double deriv[4][3] =
    {
      { -1.0/2.0,  1.0/2.0, 0.0     },
      { -1.0/2.0,  1.0/2.0, 0.0     },
      {      0.0, -1.0/5.0, 1.0/5.0 },
      {      0.0, -1.0/5.0, 1.0/5.0 }
    };

    gsl_bspline_workspace *w = gsl_bspline_alloc(2, 3);
    const size_t order = gsl_bspline_order(w);
    const size_t ncontrol = gsl_bspline_ncontrol(w);
    gsl_matrix *dB = gsl_matrix_alloc(ncontrol, order + 1);
    gsl_vector *breakpts = gsl_vector_alloc(3);
    gsl_vector_set(breakpts, 0, 0.0);
    gsl_vector_set(breakpts, 1, 2.0);
    gsl_vector_set(breakpts, 2, 7.0);
    gsl_bspline_init_augment(breakpts, w);

    for (i = 0; i < 4; ++i)  /* at each location */
      {
        /* Initialize dB with poison to ensure we overwrite it */
        gsl_matrix_set_all(dB, GSL_NAN);

        gsl_bspline_eval_deriv_basis(xloc[i], order, dB, w);
        for (j = 0; j < ncontrol; ++j)
          {
            /* check basis function 1st deriv */
            gsl_test_abs(gsl_matrix_get(dB, j, 1), deriv[i][j], GSL_DBL_EPSILON,
                         "b-spline k=%d basis #%d derivative %d at x = %f",
                         order, j, 1, xloc[i]);
          }
        for (j = 0; j < ncontrol; ++j)
          {
            /* check k order basis function has k-th deriv equal to 0 */
            gsl_test_abs(gsl_matrix_get(dB, j, order), 0.0,
                         GSL_DBL_EPSILON,
                         "b-spline k=%d basis #%d derivative %d at x = %f",
                         order, j, order,
                         xloc[i]);
          }
      }

    gsl_matrix_free(dB);
    gsl_bspline_free(w);
    gsl_vector_free(breakpts);
  }

  /* Spot check known 0th, 1st, 2nd derivative
     evaluations for a particular k = 3 case.  */
  {
    size_t i, j; /* looping */
    const double xloc[5]     =  { 0.0, 5.0, 9.0, 12.0, 15.0 };
    const double eval[5][6]  =
    {
      { 4./25.,  69./100.,   3./ 20. ,  0.    , 0.   , 0.    },
      { 0.     ,  4./21. , 143./210. ,  9./70., 0.   , 0.    },
      { 0.     ,  0.     ,   3./ 10. ,  7./10., 0.   , 0.    },
      { 0.     ,  0.     ,   0.      ,  3./4. , 1./4., 0.    },
      { 0.     ,  0.     ,   0.      ,  1./3. , 5./9., 1./9. }
    };
    const double deriv[5][6] =
    {
      { -4./25.,  3./50.,   1./ 10.,  0.    , 0.    , 0.      },
      {  0.    , -2./21.,   1./105.,  3./35., 0.    , 0.      },
      {  0.    ,  0.    ,  -1./5.  ,  1./ 5., 0.    , 0.      },
      {  0.    ,  0.    ,   0.     , -1./ 6., 1./6. , 0.      },
      {  0.    ,  0.    ,   0.     , -1./ 9., 1./27., 2./27. }
    };
    const double deriv2[5][6] =
    {
      { 2./25., -17./150.,   1.0/30.0 ,  0.0     ,  0.     , 0.     },
      { 0.    ,   1./ 42., -11.0/210.0,  1.0/35.0,  0.     , 0.     },
      { 0.    ,   0.     ,   1.0/15.0 ,-11.0/90.0,  1./18. , 0.     },
      { 0.    ,   0.     ,   0.0      ,  1.0/54.0, -7./162., 2./81. },
      { 0.    ,   0.     ,   0.0      ,  1.0/54.0, -7./162., 2./81. }
    };

    gsl_bspline_workspace *w = gsl_bspline_alloc(3, 5);
    const size_t order = gsl_bspline_order(w);
    const size_t ncontrol = gsl_bspline_ncontrol(w);
    gsl_matrix *dB = gsl_matrix_alloc(ncontrol, order + 1);
    gsl_vector *breakpts = gsl_vector_alloc(5);
    gsl_vector_set(breakpts, 0, -3.0);
    gsl_vector_set(breakpts, 1,  2.0);
    gsl_vector_set(breakpts, 2,  9.0);
    gsl_vector_set(breakpts, 3, 12.0);
    gsl_vector_set(breakpts, 4, 21.0);
    gsl_bspline_init_augment(breakpts, w);

    for (i = 0; i < 5; ++i)  /* at each location */
      {
        /* Initialize dB with poison to ensure we overwrite it */
        gsl_matrix_set_all(dB, GSL_NAN);
        gsl_bspline_eval_deriv_basis(xloc[i], order, dB, w);

        /* check basis function evaluation */
        for (j = 0; j < ncontrol; ++j)
          {
            gsl_test_abs(gsl_matrix_get(dB, j, 0), eval[i][j], GSL_DBL_EPSILON,
                         "b-spline k=%d basis #%d derivative %d at x = %f",
                         order, j, 0, xloc[i]);
          }
        /* check 1st derivative evaluation */
        for (j = 0; j < ncontrol; ++j)
          {
            gsl_test_abs(gsl_matrix_get(dB, j, 1), deriv[i][j], GSL_DBL_EPSILON,
                         "b-spline k=%d basis #%d derivative %d at x = %f",
                         order, j, 1, xloc[i]);
          }
        /* check 2nd derivative evaluation */
        for (j = 0; j < ncontrol; ++j)
          {
            gsl_test_abs(gsl_matrix_get(dB, j, 2), deriv2[i][j], GSL_DBL_EPSILON,
                         "b-spline k=%d basis #%d derivative %d at x = %f",
                         order, j, 2, xloc[i]);
          }
      }

    gsl_matrix_free(dB);
    gsl_bspline_free(w);
    gsl_vector_free(breakpts);
  }

  gsl_test(test_gram(rng_p),          "Gram Matrix");
  gsl_test(test_greville(),           "Greville");
  gsl_test(test_integ(rng_p),         "Integration");
  gsl_test(test_interp(),             "Interpolation");
  gsl_test(test_LS(rng_p),            "Least Squares");
  gsl_test(test_PLS(rng_p),           "Periodic Least Squares");

  gsl_rng_free(rng_p);

  exit(gsl_test_summary());
}
