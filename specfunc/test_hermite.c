/* specfunc/test_hermite.c
 * 
 * Copyright (C) 2011, 2012, 2013, 2014 Konrad Griessinger
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
#include <gsl/gsl_sf.h>

#include "test_sf.h"

#define WKB_TOL (1.0e+04 * TEST_SQRT_TOL0)

/*
 * Test the identities:
 *
 * Sum_{k=0}^n (n choose k) (2 y)^{n - k} H_k(x) = H_n(x + y)
 *
 * and
 *
 * Sum_{k=0}^n (n choose k) y^{n-k} He_k(x) = He_n(x + y)
 *
 * see: http://mathworld.wolfram.com/HermitePolynomial.html (Eq. 55)
 */
void
test_hermite_id1(const int n, const double x, const double y)
{
  double *a = malloc((n + 1) * sizeof(double));
  double *b = malloc((n + 1) * sizeof(double));
  double lhs, rhs;
  int k;

  a[0] = gsl_pow_int(2.0 * y, n);
  b[0] = gsl_pow_int(y, n);
  for (k = 1; k <= n; ++k)
    {
      double fac = (n - k + 1.0) / (k * y);
      a[k] = 0.5 * fac * a[k - 1];
      b[k] = fac * b[k - 1];
    }

  lhs = gsl_sf_hermite_series(n, x, a);
  rhs = gsl_sf_hermite(n, x + y);
  gsl_test_rel(lhs, rhs, TEST_TOL4, "identity1 phys n=%d x=%g y=%g", n, x, y);

  lhs = gsl_sf_hermite_prob_series(n, x, b);
  rhs = gsl_sf_hermite_prob(n, x + y);
  gsl_test_rel(lhs, rhs, TEST_TOL3, "identity1 prob n=%d x=%g y=%g", n, x, y);

  free(a);
  free(b);
}

int
test_hermite(void)
{
  gsl_sf_result r;
  
  int s = 0;
  int m, n, sa;
  double res[256];
  double x;
  const double aizero1 = -2.3381074104597670384891972524467; /* first zero of the Airy function Ai */

  /* test some known identities */
  test_hermite_id1(10, 0.75, 0.33);
  test_hermite_id1(8, -0.75, 1.20);
  test_hermite_id1(7, 2.88, -3.2);

  TEST_SF(s, gsl_sf_hermite_prob_e, (0, 0.75, &r),  1.,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_prob_e, (1, 0.75, &r),  0.75,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_prob_e, (25, 0., &r),  0.,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_prob_e, (28, 0., &r),  2.13458046676875e14,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_prob_e, (30, 0., &r), -6.190283353629375e15,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_prob_e, (25, 0.75, &r),  -1.08128685847680748265939328423e12,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_prob_e, (28, 0.75, &r),  -1.60620252094658918105511125135e14,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF_RETURN(s, gsl_sf_hermite_prob_e, (2800, 0.75, &r), GSL_EOVRFLW);
  TEST_SF_RETURN(s, gsl_sf_hermite_prob_e, (2800, 1.75, &r), GSL_EOVRFLW);

  TEST_SF(s, gsl_sf_hermite_prob_deriv_e, (225, 128, 0.75, &r),  0.,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_prob_deriv_e, (5, 128, 0.75, &r),  -3.0288278964712702882066404e112,  TEST_TOL1, GSL_SUCCESS);

  x = 0.75;
  sa = 0;
  gsl_sf_hermite_prob_array(0, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_array(0, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array(1, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[1],   +0.0,   x,   TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_array(1, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array(100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[1],   +0.0,   x, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   823.810509681701660156250, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,  1.03749254986255755872498e78, TEST_TOL1);
  gsl_test(sa, "gsl_sf_hermite_prob_array(100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array_deriv(0, 100, 0.75, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   823.810509681701660156250, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,  1.03749254986255755872498e78, TEST_TOL1);
  gsl_test(sa, "gsl_sf_hermite_prob_array_deriv(0, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array_deriv(1000, 100, 0.75, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   0.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_array_deriv(1000, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array_deriv(99, 100, 0.75, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[99], +0.0,    9.332621544394415268169923886e155, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   6.999466158295811451127442914e157, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_array_deriv(99, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array_deriv(100, 100, 0.75, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   9.332621544394415268169923886e157, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_array_deriv(100, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array_deriv(23, 100, 0.75, res);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[23], +0.0,  2.5852016738884976640000000000e22, TEST_TOL0);
  TEST_SF_VAL(sa, res[24], +0.0,  4.6533630129992957952000000000e23, TEST_TOL0);
  TEST_SF_VAL(sa, res[37], +0.0,  2.3592417210568968566591219172e37, TEST_TOL0);
  TEST_SF_VAL(sa, res[60], +0.0, -9.2573208827175536243052086845e59, TEST_TOL0);
  TEST_SF_VAL(sa, res[61], +0.0,  5.6259625429686219137261735305e59, TEST_TOL1);
  TEST_SF_VAL(sa, res[100], +0.0, 2.6503570965896336273549197e100, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_array_deriv(23, 37, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_deriv_array(100, 50, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   -3.88338863813139372375411561e31, TEST_TOL0);
  TEST_SF_VAL(sa, res[1],   +0.0,   -4.04757862431646677625108652e32, TEST_TOL1);
  TEST_SF_VAL(sa, res[10],  +0.0,    7.9614368698398116765703194e38, TEST_TOL0);
  TEST_SF_VAL(sa, res[30],  +0.0,   -9.1416928183915197188795338e54, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,    0.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_deriv_array(100, 50, 0.75)");
  s += sa;

  n = 128;
  res[0] = 1.;
  for(m=1; m<=n; m++){
    res[m] = res[m-1]/2.;
  }
  TEST_SF(s, gsl_sf_hermite_prob_series_e, (n, x, res, &r),  -4.0451066556993485405907339548e68,  TEST_TOL0, GSL_SUCCESS);

  /* phys */

  x = 0.75;

  TEST_SF(s, gsl_sf_hermite_e, (0, 0.75, &r),  1.0,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_e, (1, 0.75, &r),  1.5,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_e, (25, 0., &r),  0.,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_e, (28, 0., &r),  3.497296636753920000e18,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_e, (30, 0., &r), -2.028432049317273600e20,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_e, (25, 0.75, &r),  -9.7029819451106077507781088352e15,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_e, (28, 0.75, &r),   3.7538457078067672096408339776e18,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF_RETURN(s, gsl_sf_hermite_e, (2800, 0.75, &r), GSL_EOVRFLW);
  TEST_SF_RETURN(s, gsl_sf_hermite_e, (2800, 10.1, &r), GSL_EOVRFLW);

  TEST_SF(s, gsl_sf_hermite_deriv_e, (225, 128, 0.75, &r),  0.,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_deriv_e, (5, 128, 0.75, &r),  2.89461215568095657569833e132,  TEST_TOL0, GSL_SUCCESS);

  x = 0.75;
  sa = 0;
  gsl_sf_hermite_array(0, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_array(0, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_array(1, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[1],   +0.0,   2.*x,   TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_array(1, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_array(100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[1],   +0.0,   2.0*x, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   38740.4384765625, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,  -1.4611185395125104593177790757e93, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_array(100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_array_deriv(0, 100, 0.75, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   38740.4384765625, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,  -1.4611185395125104593177790757e93, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_array_deriv(0, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_array_deriv(1000, 100, 0.75, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   0.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_array_deriv(1000, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_array_deriv(99, 100, 0.75, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[99], +0.0,    5.915251651227242890408570128e185, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   8.872877476840864335612855192e187, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_array_deriv(99, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_array_deriv(100, 100, 0.75, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   1.183050330245448578081714026e188, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_array_deriv(100, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_array_deriv(23, 100, 0.75, res);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[23],  +0.0,   2.168624344319444261221171200e29, TEST_TOL0);
  TEST_SF_VAL(sa, res[24],  +0.0,   7.807047639549999340396216320e30, TEST_TOL0);
  TEST_SF_VAL(sa, res[37],  +0.0,   1.930387357696033719818118732e46, TEST_TOL0);
  TEST_SF_VAL(sa, res[60],  +0.0,   6.775378005383186748501182409e71, TEST_TOL0);
  TEST_SF_VAL(sa, res[61],  +0.0,  -4.451215867508936256056845902e73, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   2.957966000491202678467161e118, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_array_deriv(23, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_deriv_array(100, 50, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,  -8.26632218305863100726861832e38, TEST_TOL0);
  TEST_SF_VAL(sa, res[1],   +0.0,   2.40954750392844799126151557e40, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   1.52281030265187793605875604e49, TEST_TOL0);
  TEST_SF_VAL(sa, res[30],  +0.0,   9.52199419132990437892101915e65, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   0.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_deriv_array(100, 50, 0.75)");
  s += sa;

  n = 128;
  /* arbitrary weights */
  res[0] = 1.;
  for(m=1; m<=n; m++){
    res[m] = res[m-1]/2.;
  }
  TEST_SF(s, gsl_sf_hermite_series_e, (n, x, res, &r),  1.07772223811696567390619566842e88,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (0,  1.3, &r),  0.322651504564963746879400858624,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (0,  1.3, &r),  0.322651504564963746879400858624,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (1,  1.3, &r),  0.593187573778613235895531272243,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (1,  1.3, &r),  0.593187573778613235895531272243,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (1, -1.3, &r), -0.593187573778613235895531272243,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (1, -1.3, &r), -0.593187573778613235895531272243,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (27, 0, &r),  0.0,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (27, 0, &r),  0.0,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (28, 0, &r),  0.290371943657199641200016132937,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (28, 0, &r),  0.290371943657199641200016132937,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (28, 0.75, &r),  0.23526280808621240649319140441,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (28, 0.75, &r),  0.23526280808621240649319140441,  TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (200, 0.75, &r), -0.13725356483699291817038427801, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (200, 0.75, &r), -0.13725356483699291817038427801, TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (100028, 0.75, &r),  -0.02903467369856961147236598086,  TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (100028, 0.75, &r),  -0.02903467369856961147236598086,  TEST_TOL5, GSL_SUCCESS);

  n = 10025;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2.5);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  -0.05301278920004176459797680403,  TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (n, x, &r),  -0.05301278920004176459797680403,  TEST_TOL4, GSL_SUCCESS);

  n = 10028;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2.5);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  0.06992968509693993526829596970,  TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (n, x, &r),  0.06992968509693993526829596970,  TEST_TOL3, GSL_SUCCESS);

  n = 10025;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2.5);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  0.08049000991742150521366671021,  TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (n, x, &r),  0.08049000991742150521366671021,  TEST_TOL4, GSL_SUCCESS);

  n = 10028;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2.5);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  0.08048800667512084252723933250,  TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (n, x, &r),  0.08048800667512084252723933250,  TEST_TOL4, GSL_SUCCESS);

  n = 10025;
  x = (sqrt(2*n+1.)-2.5*(aizero1/pow(8.*n,1/6.)));
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  7.97206830806663013555068100e-6,  TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (n, x, &r),  7.97206830806663013555068100e-6,  1.0e-9, GSL_SUCCESS);

  n = 10028;
  x = (sqrt(2*n+1.)-2.5*(aizero1/pow(8.*n,1/6.)));
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  7.97188517397786729928465829e-6,  TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (n, x, &r),  7.97188517397786729928465829e-6,  1.0e-8, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_e, (10000, 60.0, &r), 0.03162606955427450540143292572, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_fast_e, (10000, 60.0, &r), 0.03162606955427450540143292572, TEST_TOL4, GSL_SUCCESS);

  x = 0.75;

  sa = 0;
  gsl_sf_hermite_func_array(100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.566979307027693616978839335983, TEST_TOL0);
  TEST_SF_VAL(sa, res[1],   +0.0,   0.601372369187597546203014795470, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.360329854170806945032958735574, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   -0.07616422890563462489003733382, TEST_TOL1);
  gsl_test(sa, "gsl_sf_hermite_func_array(100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_func_array(250, 20.0, res);
  TEST_SF_VAL(sa, res[199], +0.0,   0.254621970261340422399150964119, TEST_TOL2);
  TEST_SF_VAL(sa, res[200], +0.0,   0.288364948026205642290411878357, TEST_TOL2);
  TEST_SF_VAL(sa, res[210], +0.0,   0.266625427572822255774117774166, TEST_TOL2);
  TEST_SF_VAL(sa, res[220], +0.0,  -0.296028413764592652216845612832, TEST_TOL2);
  TEST_SF_VAL(sa, res[230], +0.0,   0.229657495510699141971182819560, TEST_TOL2);
  TEST_SF_VAL(sa, res[240], +0.0,  -0.024027870622288819185159269632, TEST_TOL2);
  TEST_SF_VAL(sa, res[250], +0.0,  -0.236101289420398152417654177998, TEST_TOL2);
  gsl_test(sa, "gsl_sf_hermite_func_array(250, 20.0)");
  s += sa;

  n = 128;
  /* arbitrary weights */
  res[0] = 1.;
  for(m=1; m<=n; m++){
    res[m] = res[m-1]/2.;
  }
  TEST_SF(s, gsl_sf_hermite_func_series_e, (n, x, res, &r),  0.81717103529960997134154552556,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_der_e, (0, 28, 0.75, &r),  0.235262808086212406493191404,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (1, 28, 0.75, &r),  1.289485094958329643927802330,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (2, 28, 0.75, &r), -13.27764473136561269145948989,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (3, 28, 0.75, &r), -72.42242083458141066943555691,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (4, 28, 0.75, &r),  753.6960554274941800190147503,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (5, 28, 0.75, &r),  4035.32788513029308540826835,   TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_func_der_e, (0, 380, 0.75, &r), -0.0400554661321992411631174,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (1, 380, 0.75, &r), -4.0417244263030600591206553,    TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (2, 380, 0.75, &r),  30.4596785269042604519780923,   TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (3, 380, 0.75, &r),  3073.4187352276349348458186556, TEST_TOL2, GSL_SUCCESS);

  {
    /* positive zeros of the probabilists' Hermite polynomial of order 17 */
    double He17z[8] = { 0.751842600703896170737870774614, 1.50988330779674075905491513417, 2.28101944025298889535537879396,
                        3.07379717532819355851658337833, 3.90006571719800990903311840097, 4.77853158962998382710540812497,
                        5.74446007865940618125547815768,6.88912243989533223256205432938 };

    n = 17;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_prob_zero_e, (n, m, &r),  He17z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  {
    /* positive zeros of the probabilists' Hermite polynomial of order 18 */
    double He18z[9] = { 0.365245755507697595916901619097, 1.09839551809150122773848360538, 1.83977992150864548966395498992,
                        2.59583368891124032910545091458, 3.37473653577809099529779309480, 4.18802023162940370448450911428,
                        5.05407268544273984538327527397, 6.00774591135959752029303858752, 7.13946484914647887560975631213 };

    n = 18;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_prob_zero_e, (n, m, &r),  He18z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  {
    /* positive zeros of the probabilists' Hermite polynomial of order 23 */
    double He23z[11] = { 0.648471153534495816722576841197, 1.29987646830397886997876116860, 1.95732755293342410739100839243,
                         2.62432363405918177067330340783, 3.30504002175296456723204112903, 4.00477532173330406712238633738,
                         4.73072419745147329426707133987, 5.49347398647179412855289367747, 6.31034985444839982842886078177,
                         7.21465943505186138595859194492, 8.29338602741735258945596770157 };

    n = 23;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_prob_zero_e, (n, m, &r),  He23z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  {
    /* positive zeros of the probabilists' Hermite polynomial of order 24 */
    double He24z[12] = { 0.317370096629452319318170455994, 0.953421922932109084904629632351, 1.59348042981642010695074168129,
                         2.24046785169175236246653790858, 2.89772864322331368932008199475, 3.56930676407356024709649151613,
                         4.26038360501990548884317727406, 4.97804137463912033462166468006, 5.73274717525120114834341822330,
                         6.54167500509863444148277523364, 7.43789066602166310850331715501, 8.50780351919525720508386233432 };

    n = 24;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_prob_zero_e, (n, m, &r),  He24z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  {
    /* positive zeros of the physicists' Hermite polynomial of order 17 */
    double H17z[8] = { 0.531633001342654731349086553718, 1.06764872574345055363045773799, 1.61292431422123133311288254454,
                       2.17350282666662081927537907149, 2.75776291570388873092640349574, 3.37893209114149408338327069289,
                       4.06194667587547430689245559698, 4.87134519367440308834927655662 };

    n = 17;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_zero_e, (n, m, &r),  H17z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  {
    /* positive zeros of the physicists' Hermite polynomial of order 18 */
    double H18z[9] = { 0.258267750519096759258116098711, 0.776682919267411661316659462284, 1.30092085838961736566626555439,
                       1.83553160426162889225383944409, 2.38629908916668600026459301424, 2.96137750553160684477863254906,
                       3.57376906848626607950067599377, 4.24811787356812646302342016090, 5.04836400887446676837203757885 };

    n = 18;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_zero_e, (n, m, &r),  H18z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  {
    /* positive zeros of the physicists' Hermite polynomial of order 23 */
    double H23z[11] = { 0.458538350068104797757887329284, 0.919151465442563765431719239593, 1.38403958568249523732634717118,
                        1.85567703767137106251504753718, 2.33701621147445578644623502174, 2.83180378712615690144806140734,
                        3.34512715994122457247439814585, 3.88447270810610186607248760288, 4.46209117374000667673186157071,
                        5.10153461047667712968749766165, 5.86430949898457256538748413474 };

    n = 23;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_zero_e, (n, m, &r),  H23z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  {
    /* positive zeros of the physicists' Hermite polynomial of order 24 */
    double H24z[12] = { 0.224414547472515585151136715527, 0.674171107037212236000245923730, 1.12676081761124507213306126773,
                        1.58425001096169414850563336202, 2.04900357366169891178708399532, 2.52388101701142697419907602333,
                        3.01254613756556482565453858421, 3.52000681303452471128987227609, 4.05366440244814950394766297923,
                        4.62566275642378726504864923776, 5.25938292766804436743072304398, 6.01592556142573971734857350899 };

    n = 24;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_zero_e, (n, m, &r),  H24z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  n = 2121;
  x = -1.*n;
  res[0] = (double) n;
  sa = 0;
  for (m=1; m<=n/2; m++) {
    gsl_sf_hermite_prob_zero_e(n, m, &r);
    if (x>=r.val) {
      sa += TEST_SF_INCONS;
      printf("sanity check failed! (gsl_sf_hermite_prob_zero)\n");
    }
    res[0] = GSL_MIN(res[0],fabs(x-r.val));
    x = r.val;
  }
  gsl_test(sa, "gsl_sf_hermite_prob_zero(n, m, r)");

  n = 2121;
  x = -1.*n;
  res[0] = (double) n;
  sa = 0;
  for (m=1; m<=n/2; m++) {
    gsl_sf_hermite_zero_e(n, m, &r);
    if (x>=r.val) {
      sa += TEST_SF_INCONS;
      printf("sanity check failed! (gsl_sf_hermite_zero)\n");
    }
    res[0] = GSL_MIN(res[0],fabs(x-r.val));
    x = r.val;
  }
  gsl_test(sa, "gsl_sf_hermite_zero(n, m, r)");

  return s;
}
