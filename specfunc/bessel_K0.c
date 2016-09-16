/* specfunc/bessel_K0.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2016 Pavel Holoborodko, Patrick Alken
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

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_bessel.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/*
 Minimax rational approximation for [0,1), peak relative error = 2.04*GSL_DBL_EPSILON.
 Source: http://www.advanpix.com/?p=3812
*/
static double k0_poly[8] = {
   1.1593151565841244842077226e-01,
   2.7898287891460317300886539e-01,
   2.5248929932161220559969776e-02,
   8.4603509072136578707676406e-04,
   1.4914719243067801775856150e-05,
   1.6271068931224552553548933e-07,
   1.2082660336282566759313543e-09,
   6.6117104672254184399933971e-12
};

static double i0_poly[7] = {
   1.0000000000000000044974165e+00,
   2.4999999999999822316775454e-01,
   2.7777777777892149148858521e-02,
   1.7361111083544590676709592e-03,
   6.9444476047072424198677755e-05,
   1.9288265756466775034067979e-06,
   3.9908220583262192851839992e-08
};

/*
 Chebyshev expansion for [1,8], peak relative error = 1.28*GSL_DBL_EPSILON. 
 Source: Pavel Holoborodko.
*/
static double ak0_data[24] = {
  -3.28737867094650101e-02,
  -4.49369057710236880e-02,
  +2.98149992004308095e-03,
  -3.03693649396187920e-04,
  +3.91085569307646836e-05,
  -5.86872422399215952e-06,
  +9.82873709937322009e-07,
  -1.78978645055651171e-07,
  +3.48332306845240957e-08,
  -7.15909210462546599e-09,
  +1.54019930048919494e-09,
  -3.44555485579194210e-10,
  +7.97356101783753023e-11,
  -1.90090968913069735e-11,
  +4.65295609304114621e-12,
  -1.16614287433470780e-12,
  +2.98554375218596891e-13,
  -7.79276979512292169e-14,
  +2.07027467168948402e-14,
  -5.58987860393825313e-15,
  +1.53202965950646914e-15,
  -4.25737536712188186e-16,
  +1.19840238501357389e-16,
  -3.41407346762502397e-17
};

static cheb_series ak0_cs = {
  ak0_data,
  23,
  -1, 1,
  10
};

/* 
 Chebyshev expansion for [8,inf), peak relative error = 1.25*GSL_DBL_EPSILON.
 Source: SLATEC/dbsk0e.f
*/
static double ak02_data[14] = {
  -.1201869826307592240E-1,
  -.9174852691025695311E-2,
  +.1444550931775005821E-3,
  -.4013614175435709729E-5,
  +.1567831810852310673E-6,
  -.7770110438521737710E-8,
  +.4611182576179717883E-9,
  -.3158592997860565771E-10,
  +.2435018039365041128E-11,
  -.2074331387398347898E-12,
  +.1925787280589917085E-13,
  -.1927554805838956104E-14,
  +.2062198029197818278E-15,
  -.2341685117579242403E-16
};

static cheb_series ak02_cs = {
  ak02_data,
  13,
  -1, 1,
  8
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K0_scaled_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 1.0) {
    const double lx = log(x);
    const double ex = exp(x);
    const double x2 = x*x;
    result->val  = ex * (gsl_poly_eval(k0_poly,8,x2)-lx*(1.0+0.25*x2*gsl_poly_eval(i0_poly,7,0.25*x2)));
    result->err  = ex * (1.6+fabs(lx)*0.6) * GSL_DBL_EPSILON;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x <= 8.0) {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak0_cs, (16.0/x-9.0)/7.0, &c);
    result->val  = (1.203125 + c.val) / sx; /* 1.203125 = 77/64 */
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak02_cs, 16.0/x-1.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = (c.err + GSL_DBL_EPSILON) / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  } 
}


int gsl_sf_bessel_K0_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 1.0) {
    const double lx = log(x);
    const double x2 = x*x;
    result->val  = gsl_poly_eval(k0_poly,8,x2)-lx*(1.0+0.25*x2*gsl_poly_eval(i0_poly,7,0.25*x2));
    result->err  = (1.6+fabs(lx)*0.6) * GSL_DBL_EPSILON;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result K0_scaled;
    int stat_K0 = gsl_sf_bessel_K0_scaled_e(x, &K0_scaled);
    int stat_e  = gsl_sf_exp_mult_err_e(-x, GSL_DBL_EPSILON*fabs(x),
                                           K0_scaled.val, K0_scaled.err,
                                           result);
    return GSL_ERROR_SELECT_2(stat_e, stat_K0);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_bessel_K0_scaled(const double x)
{
  EVAL_RESULT(gsl_sf_bessel_K0_scaled_e(x, &result));
}

double gsl_sf_bessel_K0(const double x)
{
  EVAL_RESULT(gsl_sf_bessel_K0_e(x, &result));
}

