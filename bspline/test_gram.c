/* bspline/test_gram.c
 *
 * Copyright (C) 2019, 2020, 2021 Patrick Alken
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

typedef struct
{
  size_t i;
  size_t j;
  size_t nderiv;
  gsl_bspline_workspace * w;
} numint_params;

static double
numint_func (double x, void * params)
{
  numint_params * p = (numint_params *) params;
  gsl_bspline_workspace * w = p->w;
  const size_t order = w->spline_order;
  const size_t i = p->i;
  const size_t j = p->j;
  size_t istart;

  gsl_bspline_basis_deriv(x, p->nderiv, w->dB, &istart, w);

  if (i < istart || j < istart)
    {
      return 0.0;
    }
  else if (i >= istart + order || j >= istart + order)
    {
      return 0.0;
    }
  else
    {
      double Ni = gsl_matrix_get(w->dB, i - istart, p->nderiv);
      double Nj = gsl_matrix_get(w->dB, j - istart, p->nderiv);
      return (Ni * Nj);
    }
}

/*
Calculate elements of Gram matrix using Gauss-Legendre integration on
each knot interval and compare with more efficient routine
*/

static int
test_gram_eps(const double tol, const size_t order,
              const size_t nbreak, const int interval,
              double lower, double upper, gsl_rng * rng_p)
{
  int status = GSL_SUCCESS;
  gsl_bspline_workspace * w = gsl_bspline_alloc(order, nbreak);
  const size_t ncontrol = gsl_bspline_ncontrol(w);
  gsl_vector * bkpts = gsl_vector_alloc(nbreak);
  gsl_integration_glfixed_table * gltable = gsl_integration_glfixed_table_alloc(order);
  gsl_matrix * G = gsl_matrix_alloc(ncontrol, order);
  gsl_function F;
  numint_params params;
  size_t i, j, nderiv;

  random_vector(-10.0, 10.0, bkpts, rng_p);
  gsl_sort_vector(bkpts);
  gsl_bspline_init_augment(bkpts, w);

  if (interval == 0)
    {
      lower = gsl_vector_get(w->knots, 0);
      upper = gsl_vector_get(w->knots, w->knots->size - 1);
    }

  F.function = numint_func;
  F.params = &params;

  for (nderiv = 0; nderiv < order; ++nderiv)
    {
      if (interval)
        gsl_bspline_gram_interval(lower, upper, nderiv, G, w);
      else
        gsl_bspline_gram(nderiv, G, w);

      for (j = 0; j < ncontrol; ++j)
        {
          for (i = 0; i < order && i + j < ncontrol; ++i)
            {
              double Gji = gsl_matrix_get(G, j, i);
              double Hji = 0.0;
              size_t m;

              params.i = i + j;
              params.j = j;
              params.nderiv = nderiv;
              params.w = w;

              for (m = 0; m < w->knots->size - 1; ++m)
                {
                  double a = GSL_MAX(gsl_vector_get(w->knots, m), lower);
                  double b = GSL_MIN(gsl_vector_get(w->knots, m + 1), upper);

                  if (b < a)
                    continue;

                  Hji += gsl_integration_glfixed(&F, a, b, gltable);
                }

              gsl_test_rel(Gji, Hji, tol,
                           "gram[%g,%g] random order=%zu nbreak=%zu nderiv=%zu row=%zu col=%zu",
                           lower, upper, order, nbreak, nderiv, i + j, j);
            }
        }
    }

  gsl_integration_glfixed_table_free(gltable);
  gsl_vector_free(bkpts);
  gsl_matrix_free(G);
  gsl_bspline_free(w);

  return status;
}

static int
test_gram(gsl_rng * rng_p)
{
  int status = GSL_SUCCESS;
  const double tol = 1.0e-10;
  const size_t order = 4;
  const size_t ncontrol = 7;
  gsl_bspline_workspace *w = gsl_bspline_alloc_ncontrol(order, ncontrol);
  gsl_matrix * G = gsl_matrix_alloc(ncontrol, order);
  gsl_matrix * dG = gsl_matrix_alloc(ncontrol, order);
  gsl_matrix * ddG = gsl_matrix_alloc(ncontrol, order);
  size_t i, j;

  /* these values were computed with Mathematica via numerical integration:
   * > knots = { 0,0,0,0,.25,.5,.75,1,1,1,1} 
   * > B[n_, x_] = D[BSplineBasis[{3, knots}, n, x], {x, nderiv}] 
   * > NIntegrate[B[i, x]*B[j, x], {x, 0, 1}, WorkingPrecision -> 16, PrecisionGoal -> 16]
   */
  const double G0_expected[7][4] =
  {
    {  3.571428571428571e-02,  2.187500000000000e-02,  4.613095238095238e-03,  2.976190476190476e-04, },
    {  5.535714285714286e-02,  3.906250000000000e-02,  8.630952380952381e-03,  7.440476190476190e-05 },
    {  8.169642857142857e-02,  5.615079365079365e-02,  5.902777777777778e-03,  7.440476190476190e-05 },
    {  1.198412698412698e-01,  5.615079365079365e-02,  8.630952380952381e-03,  2.976190476190476e-04 },
    {  8.169642857142857e-02,  3.906250000000000e-02,  4.613095238095238e-03,  0.000000000000000e+00 },
    {  5.535714285714286e-02,  2.187500000000000e-02,  0.000000000000000e+00,  0.000000000000000e+00 },
    {  3.571428571428571e-02,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00 }
  };

  const double dG0_expected[7][4] =
  {
    {  7.200000000000000e+00, -5.100000000000000e+00, -1.900000000000000e+00, -2.000000000000000e-01 },
    {  6.000000000000000e+00,  1.500000000000000e-01, -1.000000000000000e+00, -5.000000000000000e-02 },
    {  2.700000000000000e+00, -1.333333333333333e-01, -7.666666666666667e-01, -5.000000000000000e-02 },
    {  2.666666666666667e+00, -1.333333333333333e-01, -1.000000000000000e+00, -2.000000000000000e-01 },
    {  2.700000000000000e+00,  1.500000000000000e-01, -1.900000000000000e+00,  0.000000000000000e+00 },
    {  6.000000000000000e+00, -5.100000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00 },
    {  7.200000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00 }
  };

  const double ddG0_expected[7][4] =
  {
    {  7.680000000000000e+02, -1.056000000000000e+03,  2.240000000000000e+02,  6.400000000000000e+01 },
    {  1.536000000000000e+03, -4.320000000000000e+02, -6.400000000000000e+01,  1.600000000000000e+01 },
    {  2.880000000000000e+02, -8.533333333333333e+01, -1.066666666666667e+01,  1.600000000000000e+01 },
    {  1.706666666666667e+02, -8.533333333333333e+01, -6.400000000000000e+01,  6.400000000000000e+01 },
    {  2.880000000000000e+02, -4.320000000000000e+02,  2.240000000000000e+02,  0.000000000000000e+00 },
    {  1.536000000000000e+03, -1.056000000000000e+03,  0.000000000000000e+00,  0.000000000000000e+00 },
    {  7.680000000000000e+02,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00 }
  };

  gsl_bspline_init_uniform(0.0, 1.0, w);
  gsl_bspline_gram(0, G, w);
  gsl_bspline_gram(1, dG, w);
  gsl_bspline_gram(2, ddG, w);

  for (j = 0; j < ncontrol; ++j)
    {
      for (i = 0; i < order && i + j < ncontrol; ++i)
        {
          double Gij = gsl_matrix_get(G, j, i);
          double dGij = gsl_matrix_get(dG, j, i);
          double ddGij = gsl_matrix_get(ddG, j, i);

          gsl_test_rel(Gij, G0_expected[j][i], tol,
                       "b-spline k=%zu Gram case 0 nderiv=0 row=%zu col=%zu",
                       order, i + j, j);
          gsl_test_rel(dGij, dG0_expected[j][i], tol,
                       "b-spline k=%zu Gram case 0 nderiv=1 row=%zu col=%zu",
                       order, i + j, j);
          gsl_test_rel(ddGij, ddG0_expected[j][i], tol,
                       "b-spline k=%zu Gram case 0 nderiv=2 row=%zu col=%zu",
                       order, i + j, j);
        }
    }

  gsl_matrix_free(G);
  gsl_matrix_free(dG);
  gsl_matrix_free(ddG);
  gsl_bspline_free(w);

  /* generate random knots and test against brute force integration */
  test_gram_eps(tol, 3, 5, 0, 0.0, 0.0, rng_p);
  test_gram_eps(tol, 3, 6, 0, 0.0, 0.0, rng_p);
  test_gram_eps(tol, 4, 5, 0, 0.0, 0.0, rng_p);
  test_gram_eps(tol, 4, 10, 0, 0.0, 0.0, rng_p);
  test_gram_eps(tol, 5, 7, 0, 0.0, 0.0, rng_p);
  test_gram_eps(tol, 5, 8, 0, 0.0, 0.0, rng_p);

  test_gram_eps(tol, 3, 5, 1, 1.0, 3.0, rng_p);
  test_gram_eps(tol, 3, 6, 1, -2.0, 0.0, rng_p);
  test_gram_eps(tol, 4, 5, 1, -10.0, 9.0, rng_p);
  test_gram_eps(tol, 4, 10, 1, -1.0, 4.0, rng_p);
  test_gram_eps(tol, 5, 7, 1, 0.1, 0.2, rng_p);
  test_gram_eps(tol, 5, 8, 1, 1.0, 2.0, rng_p);

  return status;
}
