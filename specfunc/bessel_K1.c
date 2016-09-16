/* specfunc/bessel_K1.c
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
 Minimax rational approximation for [0,1), peak relative error = 1.83*GSL_DBL_EPSILON.
 Source: http://www.advanpix.com/?p=3987
*/
static double k1_poly[9] = {
  -3.0796575782920622440538935e-01,
  -8.5370719728650778045782736e-02,
  -4.6421827664715603298154971e-03,
  -1.1253607036630425931072996e-04,
  -1.5592887702110907110292728e-06,
  -1.4030163679125934402498239e-08,
  -8.8718998640336832196558868e-11,
  -4.1614323580221539328960335e-13,
  -1.5261293392975541707230366e-15
};

static double i1_poly[7] = {
  8.3333333333333325191635191e-02,
  6.9444444444467956461838830e-03,
  3.4722222211230452695165215e-04,
  1.1574075952009842696580084e-05,
  2.7555870002088181016676934e-07,
  4.9724386164128529514040614e-09
};

/*
 Chebyshev expansion for [1,8], peak relative error = 1.28*GSL_DBL_EPSILON. 
 Source: Pavel Holoborodko.
*/
static double ak1_data[25] = {
  +2.07996868001418246e-01,
  +1.62581565017881476e-01,
  -5.87070423518863640e-03,
  +4.95021520115789501e-04,
  -5.78958347598556986e-05,
  +8.18614610209334726e-06,
  -1.31604832009487277e-06,
  +2.32546031520101213e-07,
  -4.42206518311557987e-08,
  +8.92163994883100361e-09,
  -1.89046270526983427e-09,
  +4.17568808108504702e-10,
  -9.55912361791375794e-11,
  +2.25769353153867758e-11,
  -5.48128000211158482e-12,
  +1.36386122546441926e-12,
  -3.46936690565986409e-13,
  +9.00354564415705942e-14,
  -2.37950577776254432e-14,
  +6.39447503964025336e-15,
  -1.74498363492322044e-15,
  +4.82994547989290473e-16,
  -1.35460927805445606e-16,
  +3.84604274446777234e-17,
  -1.10456856122581316e-17
};

static cheb_series ak1_cs = {
  ak1_data,
  24,
  -1, 1,
  9
};

/* 
 Chebyshev expansion for [8,inf), peak relative error = 1.25*GSL_DBL_EPSILON.
 Source: SLATEC/dbsk1e.f
*/
static double ak12_data[14] = {
  +.637930834373900104E-1,
  +.283288781304972094E-1,
  -.247537067390525035E-3,
  +.577197245160724882E-5,
  -.206893921953654830E-6,
  +.973998344138180418E-8,
  -.558533614038062498E-9,
  +.373299663404618524E-10,
  -.282505196102322545E-11,
  +.237201900248414417E-12,
  -.217667738799175398E-13,
  +.215791416161603245E-14,
  -.229019693071826928E-15,
  +.258288572982327496E-16
};

static cheb_series ak12_cs = {
  ak12_data,
  13,
  -1, 1,
  7
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K1_scaled_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*GSL_DBL_MIN) {
    OVERFLOW_ERROR(result);
  }
  else if(x < 1.0) {
    const double lx = log(x);
    const double ex = exp(x);
    const double x2 = x*x;
    const double t  = 0.25*x2;    
    const double i1 = 0.5 * x * (1.0 + t * (0.5 + t * gsl_poly_eval(i1_poly,6,t)));
    result->val  = ex * (x2 * gsl_poly_eval(k1_poly,9,x2) + x * lx * i1 + 1) / x;
    result->err  = ex * (1.6+fabs(lx)*0.6) * GSL_DBL_EPSILON;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x <= 8.0) {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak1_cs, (16.0/x-9.0)/7.0, &c);
    result->val  = (1.375 + c.val) / sx; /* 1.375 = 11/8 */
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak12_cs, 16.0/x-1.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_K1_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 2.0*GSL_DBL_MIN) {
    OVERFLOW_ERROR(result);
  }
  else if(x < 1.0) {
    const double lx = log(x);
    const double x2 = x*x;
    const double t  = 0.25*x2;    
    const double i1 = 0.5 * x * (1.0 + t * (0.5 + t * gsl_poly_eval(i1_poly,6,t)));
    result->val  = (x2 * gsl_poly_eval(k1_poly,9,x2) + x * lx * i1 + 1) / x;
    result->err  = (1.6+fabs(lx)*0.6) * GSL_DBL_EPSILON;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;      
  }
  else {
    gsl_sf_result K1_scaled;
    int stat_K1 = gsl_sf_bessel_K1_scaled_e(x, &K1_scaled);
    int stat_e  = gsl_sf_exp_mult_err_e(-x, 0.0,
                                           K1_scaled.val, K1_scaled.err,
                                           result);
    result->err = fabs(result->val) * (GSL_DBL_EPSILON*fabs(x) + K1_scaled.err/K1_scaled.val);
    return GSL_ERROR_SELECT_2(stat_e, stat_K1);
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_bessel_K1_scaled(const double x)
{
  EVAL_RESULT(gsl_sf_bessel_K1_scaled_e(x, &result));
}

double gsl_sf_bessel_K1(const double x)
{
  EVAL_RESULT(gsl_sf_bessel_K1_e(x, &result));
}
