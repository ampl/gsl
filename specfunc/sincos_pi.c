/* specfunc/sincos_pi.c
 * 
 * Copyright (C) 2017 Gerard Jungman, Konrad Griessinger (konradg@gmx.net)
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

/* routines for computing sin(pi*x) and cos(pi*x), respectively, with argument reduction */

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_sincos_pi.h>

/* Any double precision number bigger than this is automatically an even integer. */
#define TWOBIG (2.0 / GSL_DBL_EPSILON)

/* routine computing sin(pi*x) valid for |x| <= 0.25 using a Taylor expansion around the origin and otherwise a rational approximation from the reference below. Spot-checked to give around 2e-16 relative accuracy. */
/* I. Koren and O. Zinaty. Evaluating elementary functions in a numerical
coprocessor based on rational approximations. IEEE Transactions on
Computers, Vol.39, No.8, August 1990, pp 1030-1037. */
/*
static int
sin_pi_koren(const double x, gsl_sf_result *result)
{
  result->val = 0.0;
  result->err = 0.0;
  if (16.0*fabs(x) < 1.0) {
    const double y = M_PI * x;
    const double a = y*y;
    result->val = y*(1.0 - a*(1.0 - a*(1.0 - a*(1.0 - a*(1.0 - a/110.0)/72.0)/42.0)/20.0)/6.0);
  }
  else {
    const double a0 = 1805490264.690988571178600370234394843221;
    const double a1 = -164384678.227499837726129612587952660511;
    const double a2 =    3664210.647581261810227924465160827365;
    const double a3 =     -28904.140246461781357223741935980097;
    const double a4 =         76.568981088717405810132543523682;
    const double b0 = 2298821602.638922662086487520330827251172;
    const double b1 =   27037050.118894436776624866648235591988;
    const double b2 =     155791.388546947693206469423979505671;
    const double b3 =        540.567501261284024767779280700089;
    const double t = 16.0*x*x;
    result->val = 4.0*x*(((( a4*t + a3 )*t + a2 )*t + a1 )*t + a0)/(((( t + b3 )*t + b2 )*t + b1 )*t + b0);
  }
  
  result->err = GSL_DBL_EPSILON*fabs(result->val);
  
  return GSL_SUCCESS;
}
*/

/* routine computing cos(pi*x) valid for |x| <= 0.25 using a Taylor expansion around the origin and otherwise a rational approximation from the reference below. Spot-checked to give around 2e-16 relative accuracy. */
/* I. Koren and O. Zinaty. Evaluating elementary functions in a numerical
coprocessor based on rational approximations. IEEE Transactions on
Computers, Vol.39, No.8, August 1990, pp 1030-1037. */
/*
static int
cos_pi_koren(const double x, gsl_sf_result *result)
{
  result->val = 0.0;
  result->err = 0.0;
  if (20.0*fabs(x) < 1.0) {
    const double y = M_PI * x;
    const double a = y*y;
    result->val = 1.0 - 0.5*a*(1.0 - a*(1.0 - a*(1.0 - a*(1.0 - a/90.0)/56.0)/30.0)/12.0);
  }
  else {
    const double a0 = 1090157078.174871420428849017262549038606;
    const double a1 = -321324810.993150712401352959397648541681;
    const double a2 =   12787876.849523878944051885325593878177;
    const double a3 =    -150026.206045948110568310887166405972;
    const double a4 =        538.333564203182661664319151379451;
    const double b0 = 1090157078.174871420428867295670039506886;
    const double b1 =   14907035.776643879767410969509628406502;
    const double b2 =     101855.811943661368302608146695082218;
    const double b3 =        429.772865107391823245671264489311;
    const double t = 16.0*x*x;
    result->val = (((( a4*t + a3 )*t + a2 )*t + a1 )*t + a0)/(((( t + b3 )*t + b2 )*t + b1 )*t + b0);
  }
  
  result->err = GSL_DBL_EPSILON*fabs(result->val);
  
  return GSL_SUCCESS;
}
*/

/* routine computing sin(pi*x) using a Taylor expansion around the origin and otherwise the library function. */
static int
sin_pi_taylor(const double x, gsl_sf_result *result)
{
  result->val = 0.0;
  result->err = 0.0;
  if (16.0*fabs(x) < 1.0) {
    const double y = M_PI * x;
    const double a = y*y;
    result->val = y*(1.0 - a*(1.0 - a*(1.0 - a*(1.0 - a*(1.0 - a/110.0)/72.0)/42.0)/20.0)/6.0);
  }
  else {
    result->val = sin(M_PI*x);
  }
  
  result->err = GSL_DBL_EPSILON*fabs(result->val);
  
  return GSL_SUCCESS;
}

/* routine computing sin(pi*x) using a Taylor expansion around the origin and otherwise the library function. */
static int
cos_pi_taylor(const double x, gsl_sf_result *result)
{
  result->val = 0.0;
  result->err = 0.0;
  if (20.0*fabs(x) < 1.0) {
    const double y = M_PI * x;
    const double a = y*y;
    result->val = 1.0 - 0.5*a*(1.0 - a*(1.0 - a*(1.0 - a*(1.0 - a/90.0)/56.0)/30.0)/12.0);
  }
  else {
    result->val = cos(M_PI*x);
  }
  
  result->err = GSL_DBL_EPSILON*fabs(result->val);
  
  return GSL_SUCCESS;
}

int
gsl_sf_sin_pi_e(const double x, gsl_sf_result *result)
{
  double intx = 0.0, fracx = 0.0;
  long q;
  int sign = 1, status;

  result->val = 0.0;
  result->err = 0.0;
  fracx = modf(x,&intx);
  if (fracx == 0.0) return GSL_SUCCESS;
  if(fabs(intx) >= TWOBIG) return GSL_SUCCESS; /* to be sure. Actually should be covered by the line above */

  q = ( ( (intx >= LONG_MIN) && (intx <= LONG_MAX) ) ? intx : fmod(intx, 2.0) );
  sign = ( q % 2 ? -1 : 1 );

  /* int sign = 1 - 2*((int)round(fmod(fabs(intx),2.0))); */
  if (fabs(fracx) == 0.5) { /* probably unnecessary */
    if (fracx < 0.0) sign = -sign;
    result->val = ( sign != 1 ? -1.0 : 1.0 );
    return GSL_SUCCESS;
  }
  if (fabs(fracx) > 0.5) {
    sign = -sign;
    fracx = ( fracx > 0.0 ? fracx-1.0 : fracx+1.0 );
  }

  status = 0;
  if (fracx > 0.25) {
    status = cos_pi_taylor((fracx-0.5), result);
  }
  else if (fracx < -0.25) {
    status = cos_pi_taylor((fracx+0.5), result);
    sign = -sign;
  }
  else {
    status = sin_pi_taylor(fracx, result);
  }
  if (sign != 1) result->val = -result->val;
  return status;
}

int
gsl_sf_cos_pi_e(const double x, gsl_sf_result *result)
{
  double intx = 0.0, fracx = 0.0;
  long q;
  int sign = 1, status;

  result->val = 0.0;
  result->err = 0.0;
  fracx = modf(x,&intx);
  if (fabs(fracx) == 0.5) return GSL_SUCCESS;
  
  if(fabs(intx) >= TWOBIG) {
    result->val = 1.0;
    return GSL_SUCCESS;
  }

  q = ( ( (intx >= LONG_MIN) && (intx <= LONG_MAX) ) ? intx : fmod(intx, 2.0) );
  sign = ( q % 2 ? -1 : 1 );

  /* int sign = 1 - 2*((int)round(fmod(fabs(intx),2.0))); */
  if (fracx == 0.0) { /* probably unnecessary */
    result->val = ( sign != 1 ? -1.0 : 1.0 );
    return GSL_SUCCESS;
  }
  if (fabs(fracx) > 0.5) {
    sign = -sign;
    fracx = ( fracx > 0.0 ? fracx-1.0 : fracx+1.0 );
  }

  status = 0;
  if (fracx > 0.25) {
    status = sin_pi_taylor((fracx-0.5), result);
    sign = -sign;
  }
  else if (fracx < -0.25) {
    status = sin_pi_taylor((fracx+0.5), result);
  }
  else {
    status = cos_pi_taylor(fracx, result);
  }
  if (sign != 1) result->val = -result->val;
  return status;
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double
gsl_sf_sin_pi(const double x)
{
  EVAL_RESULT(gsl_sf_sin_pi_e(x, &result));
}

double
gsl_sf_cos_pi(const double x)
{
  EVAL_RESULT(gsl_sf_cos_pi_e(x, &result));
}
