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

  lhs = gsl_sf_hermite_phys_series(n, x, a);
  rhs = gsl_sf_hermite_phys(n, x + y);
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

  x = 0.75;

  TEST_SF(s, gsl_sf_hermite_prob_e, (0, 0.75, &r),  1.,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_prob_e, (1, 0.75, &r),  x,  TEST_TOL0, GSL_SUCCESS);

  n = 25;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, 0., &r),  0.,  TEST_TOL0, GSL_SUCCESS);
  n = 28;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, 0., &r),  213458046676875,  TEST_TOL0, GSL_SUCCESS);
  n = 25;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, 0.75, &r),  -1.08128685847680748265939328423e12,  TEST_TOL0, GSL_SUCCESS);
  n = 28;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, 0.75, &r),  -1.60620252094658918105511125135e14,  TEST_TOL0, GSL_SUCCESS);

#if 0
  n = 10025;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2)*M_SQRT2;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, x, &r),  -3.3090527852387782540121569578e18961,  TEST_TOL0, GSL_SUCCESS);

  n = 10028;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2)*M_SQRT2;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, x, &r),  -7.515478445930242044360704363e18967,  TEST_TOL0, GSL_SUCCESS);

  n = 10025;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2)*M_SQRT2;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, x, &r),  4.1369269649092456235914193753e22243,  TEST_TOL0, GSL_SUCCESS);

  n = 10028;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2)*M_SQRT2;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, x, &r),  8.363694992558646923734666303e22250,  TEST_TOL0, GSL_SUCCESS);

  n = 10025;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)))*M_SQRT2;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, x, &r),  7.398863979737363164340057757e22273,  TEST_TOL0, GSL_SUCCESS);

  n = 10028;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)))*M_SQRT2;
  TEST_SF(s, gsl_sf_hermite_prob_e, (n, x, &r),  1.507131397474022356488976968e22281,  TEST_TOL0, GSL_SUCCESS);
#endif

  x = 0.75;

  n = 128;
  m = 225;
  TEST_SF(s, gsl_sf_hermite_prob_der_e, (m, n, x, &r),  0.,  TEST_TOL0, GSL_SUCCESS);

  n = 128;
  m = 5;
  TEST_SF(s, gsl_sf_hermite_prob_der_e, (m, n, x, &r),  -3.0288278964712702882066404e112,  TEST_TOL1, GSL_SUCCESS);

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
  TEST_SF_VAL(sa, res[10],  +0.0,   823.810509681701660156250, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,  1.03749254986255755872498e78, TEST_TOL1);
  gsl_test(sa, "gsl_sf_hermite_prob_array(100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array_der(0, 100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   823.810509681701660156250, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,  1.03749254986255755872498e78, TEST_TOL1);
  gsl_test(sa, "gsl_sf_hermite_prob_array_der(0, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array_der(1000, 100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   0.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_array_der(1000, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_array_der(23, 100, x, res);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[37], +0.0,  2.3592417210568968566591219172e37, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0, 2.6503570965896336273549197e100, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_array(23, 37, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_prob_der_array(100, 50, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   -3.88338863813139372375411561e31, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   7.9614368698398116765703194e38, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   0.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_prob_der_array(100, 50, 0.75)");
  s += sa;

  n = 128;
  res[0] = 1.;
  for(m=1; m<=n; m++){
    res[m] = res[m-1]/2.;
  }
  TEST_SF(s, gsl_sf_hermite_prob_series_e, (n, x, res, &r),  -4.0451066556993485405907339548e68,  TEST_TOL0, GSL_SUCCESS);

  /* phys */

  x = 0.75;

  TEST_SF(s, gsl_sf_hermite_phys_e, (0, 0.75, &r),  1.,  TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hermite_phys_e, (1, 0.75, &r),  2.*x,  TEST_TOL0, GSL_SUCCESS);

  n = 25;
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, 0., &r),  0.,  TEST_TOL0, GSL_SUCCESS);
  n = 28;
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, 0., &r),  3497296636753920000,  TEST_TOL0, GSL_SUCCESS);
  n = 25;
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, 0.75, &r),  -9.7029819451106077507781088352e15,  TEST_TOL0, GSL_SUCCESS);
  n = 28;
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, 0.75, &r),  3.7538457078067672096408339776e18,  TEST_TOL0, GSL_SUCCESS);

#if 0
  n = 10025;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2);
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, x, &r),  -2.7074282783315424535693575770e20470,  TEST_TOL0, GSL_SUCCESS);

  n = 10028;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2);
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, x, &r),  -1.7392214893577690864150561850e20477,  TEST_TOL0, GSL_SUCCESS);

  n = 10025;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2);
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, x, &r),  3.3847852473526979744366379542e23752,  TEST_TOL0, GSL_SUCCESS);

  n = 10028;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2);
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, x, &r),  1.9355145738418079256435712027e23760,  TEST_TOL0, GSL_SUCCESS);

  n = 10025;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)));
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, x, &r),  6.053663953512337293393128307e23782,  TEST_TOL0, GSL_SUCCESS);

  n = 10028;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)));
  TEST_SF(s, gsl_sf_hermite_phys_e, (n, x, &r),  3.487782358276961096026268141e23790,  TEST_TOL0, GSL_SUCCESS);
#endif

  x = 0.75;

  n = 128;
  m = 225;
  TEST_SF(s, gsl_sf_hermite_phys_der_e, (m, n, x, &r),  0.,  TEST_TOL0, GSL_SUCCESS);

  n = 128;
  m = 5;
  TEST_SF(s, gsl_sf_hermite_phys_der_e, (m, n, x, &r),  2.89461215568095657569833e132,  TEST_TOL0, GSL_SUCCESS);

  sa = 0;
  gsl_sf_hermite_phys_array(0, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_phys_array(0, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_phys_array(1, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[1],   +0.0,   2.*x,   TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_phys_array(1, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_phys_array(100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   38740.4384765625, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,  -1.4611185395125104593177790757e93, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_phys_array(100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_phys_array_der(0, 100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   1.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   38740.4384765625, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,  -1.4611185395125104593177790757e93, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_phys_array_der(0, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_phys_array_der(1000, 100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   0.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_phys_array_der(1000, 100, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_phys_array_der(23, 100, x, res);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.0, TEST_TOL0);
  TEST_SF_VAL(sa, res[37],  +0.0,   1.930387357696033719818118732e46, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   2.957966000491202678467161e118, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_phys_array(23, 37, 0.75)");
  s += sa;

  sa = 0;
  gsl_sf_hermite_phys_der_array(100, 50, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   -8.2663221830586310072686183e38, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   1.52281030265187793605875604e49, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   0.0, TEST_TOL0);
  gsl_test(sa, "gsl_sf_hermite_phys_der_array(100, 50, 0.75)");
  s += sa;

  n = 128;
  /* arbitrary weights */
  res[0] = 1.;
  for(m=1; m<=n; m++){
    res[m] = res[m-1]/2.;
  }
  TEST_SF(s, gsl_sf_hermite_phys_series_e, (n, x, res, &r),  1.07772223811696567390619566842e88,  TEST_TOL0, GSL_SUCCESS);

  x = 0.75;
  n = 28;
  TEST_SF(s, gsl_sf_hermite_func_e, (n, 0, &r),  0.290371943657199641200016132937,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  0.23526280808621240649319140441,  TEST_TOL0, GSL_SUCCESS);
  n = 100028;
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  -0.02903467369856961147236598086,  TEST_TOL4, GSL_SUCCESS);

  n = 10025;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2.5);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  -0.05301278920004176459797680403,  TEST_TOL4, GSL_SUCCESS);

  n = 10028;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2.5);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  0.06992968509693993526829596970,  TEST_TOL3, GSL_SUCCESS);

  n = 10025;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2.5);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  0.08049000991742150521366671021,  TEST_TOL4, GSL_SUCCESS);

  n = 10028;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2.5);
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  0.08048800667512084252723933250,  TEST_TOL4, GSL_SUCCESS);

  n = 10025;
  x = (sqrt(2*n+1.)-2.5*(aizero1/pow(8.*n,1/6.)));
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  7.97206830806663013555068100e-6,  TEST_TOL4, GSL_SUCCESS);

  n = 10028;
  x = (sqrt(2*n+1.)-2.5*(aizero1/pow(8.*n,1/6.)));
  TEST_SF(s, gsl_sf_hermite_func_e, (n, x, &r),  7.97188517397786729928465829e-6,  TEST_TOL4, GSL_SUCCESS);

  x = 0.75;

  sa = 0;
  gsl_sf_hermite_func_array(100, x, res);
  TEST_SF_VAL(sa, res[0],   +0.0,   0.566979307027693616978839335983, TEST_TOL0);
  TEST_SF_VAL(sa, res[10],  +0.0,   0.360329854170806945032958735574, TEST_TOL0);
  TEST_SF_VAL(sa, res[100], +0.0,   -0.07616422890563462489003733382, TEST_TOL1);
  gsl_test(sa, "gsl_sf_hermite_func_array(100, 0.75)");
  s += sa;

  n = 128;
  /* arbitrary weights */
  res[0] = 1.;
  for(m=1; m<=n; m++){
    res[m] = res[m-1]/2.;
  }
  TEST_SF(s, gsl_sf_hermite_func_series_e, (n, x, res, &r),  0.81717103529960997134154552556,  TEST_TOL0, GSL_SUCCESS);

  m = 5;
  n = 28;
  TEST_SF(s, gsl_sf_hermite_func_der_e, (0, n, x, &r),  0.235262808086212406493191404,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hermite_func_der_e, (m, n, x, &r),  4035.32788513029308540826835,  TEST_TOL1, GSL_SUCCESS);

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
      TEST_SF(s, gsl_sf_hermite_phys_zero_e, (n, m, &r),  H17z[m-1],  TEST_TOL0, GSL_SUCCESS);
    }
  }

  {
    /* positive zeros of the physicists' Hermite polynomial of order 18 */
    double H18z[9] = { 0.258267750519096759258116098711, 0.776682919267411661316659462284, 1.30092085838961736566626555439,
                       1.83553160426162889225383944409, 2.38629908916668600026459301424, 2.96137750553160684477863254906,
                       3.57376906848626607950067599377, 4.24811787356812646302342016090, 5.04836400887446676837203757885 };

    n = 18;
    for (m=1; m<=n/2; m++) {
      TEST_SF(s, gsl_sf_hermite_phys_zero_e, (n, m, &r),  H18z[m-1],  TEST_TOL0, GSL_SUCCESS);
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
      TEST_SF(s, gsl_sf_hermite_phys_zero_e, (n, m, &r),  H23z[m-1],  TEST_TOL0, GSL_SUCCESS);
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
      TEST_SF(s, gsl_sf_hermite_phys_zero_e, (n, m, &r),  H24z[m-1],  TEST_TOL0, GSL_SUCCESS);
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
    gsl_sf_hermite_phys_zero_e(n, m, &r);
    if (x>=r.val) {
      sa += TEST_SF_INCONS;
      printf("sanity check failed! (gsl_sf_hermite_phys_zero)\n");
    }
    res[0] = GSL_MIN(res[0],fabs(x-r.val));
    x = r.val;
  }
  gsl_test(sa, "gsl_sf_hermite_phys_zero(n, m, r)");

  return s;
}
