/* bspline/test_interp.c
 *
 * Copyright (C) 2019, 2020, 2021, 2022 Patrick Alken
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

static void
test_interp_eps (const size_t n, const size_t order, const double tol,
                 const double * data_x, const double * data_y, const char * desc)
{
  gsl_vector_const_view xv = gsl_vector_const_view_array(data_x, n);
  gsl_vector_const_view yv = gsl_vector_const_view_array(data_y, n);
  gsl_bspline_workspace * work = gsl_bspline_alloc_ncontrol(order, n);
  gsl_matrix * XB = gsl_matrix_alloc(n, 3*(order-1) + 1);
  gsl_vector * c = gsl_vector_alloc(n);
  gsl_vector_uint * ipiv = gsl_vector_uint_alloc(n);
  size_t i;

  gsl_bspline_init_interp(&xv.vector, work);

  /* test knots satisfy Schoenberg-Whitney conditions */
  for (i = 0; i < n; ++i)
    {
      double ti = gsl_vector_get(work->knots, i);
      double tipk = gsl_vector_get(work->knots, i + order);
      double xi = data_x[i];
      int s;
      
      if (i == 0 || i == n - 1)
        s = (xi < ti || xi > tipk);
      else
        s = (xi <= ti || xi >= tipk);

      gsl_test(s, "shoenberg-whitney i=%zu xi=%g [%g,%g]", i, xi, ti, tipk);
    }

  gsl_bspline_col_interp(&xv.vector, XB, work);

  /* solve linear system with banded LU */
  gsl_linalg_LU_band_decomp(n, order - 1, order - 1, XB, ipiv);
  gsl_linalg_LU_band_solve(order - 1, order - 1, XB, ipiv, &yv.vector, c);

  for (i = 0; i < n; ++i)
    {
      double result;

      gsl_bspline_calc(data_x[i], c, &result, work);
      gsl_test_rel(result, data_y[i], tol, "%s order=%zu i=%zu", desc, order, i);
    }
  
  gsl_bspline_free(work);
  gsl_matrix_free(XB);
  gsl_vector_free(c);
  gsl_vector_uint_free(ipiv);
}

static void
test_interp_hermite_eps (const size_t n, const size_t order, const double tol,
                         const double * data_x, const double * data_y, const double * data_dy,
                         const char * desc)
{
  const size_t nderiv = 1;
  gsl_vector_const_view xv = gsl_vector_const_view_array(data_x, n);
  gsl_vector_const_view yv = gsl_vector_const_view_array(data_y, n);
  gsl_vector_const_view dyv = gsl_vector_const_view_array(data_dy, n);
  const size_t ncontrol = (nderiv + 1) * (n + 2) - order;
  gsl_bspline_workspace * work = gsl_bspline_alloc_ncontrol(order, ncontrol);
  gsl_vector * c = gsl_vector_alloc(ncontrol);
  size_t i;

  gsl_bspline_init_hermite(nderiv, &xv.vector, work);
  gsl_bspline_interp_chermite(&xv.vector, &yv.vector, &dyv.vector, c, work);

  for (i = 0; i < n; ++i)
    {
      double result, dresult;

      gsl_bspline_calc(data_x[i], c, &result, work);
      gsl_bspline_calc_deriv(data_x[i], c, nderiv, &dresult, work);

      gsl_test_rel(result, data_y[i], tol, "%s order=%zu i=%zu", desc, order, i);
      gsl_test_rel(dresult, data_dy[i], tol, "%s deriv order=%zu i=%zu", desc, order, i);
    }
  
  gsl_bspline_free(work);
  gsl_vector_free(c);
}

static int
test_interp (void)
{
  int status = GSL_SUCCESS;

  {
    double data_x[6] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
    double data_y[6] = { 1.0, 0.961538461538461, 0.862068965517241, 
                         0.735294117647059, 0.609756097560976, 0.500000000000000 };
    double data_dy[6] = {  0.000000000000e+00, -3.698224852071e-01, -5.945303210464e-01,
                          -6.487889273356e-01, -5.948839976205e-01, -5.000000000000e-01 };

    test_interp_eps(6, 2, GSL_DBL_EPSILON, data_x, data_y, "1/(1+x^2) interpolation");
    test_interp_eps(6, 3, GSL_DBL_EPSILON, data_x, data_y, "1/(1+x^2) interpolation");
    test_interp_eps(6, 4, 10.0 * GSL_DBL_EPSILON, data_x, data_y, "1/(1+x^2) interpolation");

    test_interp_hermite_eps(6, 4, 1.0e2 * GSL_DBL_EPSILON, data_x, data_y, data_dy, "1/(1+x^2) Hermite interpolation");
  }

  {
    double data_x[7] = {   -1.2139767065644265, -0.792590494453907, -0.250954683125019,
                            0.665867809951305,   0.735655088722706,  0.827622053027153,
                            1.426592227816582 };
    double data_y[7] = {   -0.00453877449035645,  0.49763182550668716, 0.17805472016334534,
                            0.40514493733644485, -0.21595209836959839, 0.47405586764216423,
                            0.46561462432146072 };

    test_interp_eps(7, 2, GSL_DBL_EPSILON, data_x, data_y, "random interpolation");
    test_interp_eps(7, 3, 1.0e1 * GSL_DBL_EPSILON, data_x, data_y, "random interpolation");
    test_interp_eps(7, 4, 1.0e2 * GSL_DBL_EPSILON, data_x, data_y, "random interpolation");
    test_interp_eps(7, 5, 1.0e3 * GSL_DBL_EPSILON, data_x, data_y, "random interpolation");
  }

  return status;
}
