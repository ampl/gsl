/* specfunc/hermite.c
 * 
 * Copyright (C) 2011, 2012, 2013, 2014 Konrad Griessinger
 * (konradg(at)gmx.net)
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

/*----------------------------------------------------------------------*
 * "The purpose of computing is insight, not numbers." - R.W. Hamming   *
 * Hermite polynomials, Hermite functions                               *
 * and their respective arbitrary derivatives                           *
 *----------------------------------------------------------------------*/

/* TODO:
 * - array functions for derivatives of Hermite functions
 * - asymptotic approximation for derivatives of Hermite functions
 * - refine existing asymptotic approximations, especially around x=sqrt(2*n+1) or x=sqrt(2*n+1)*sqrt(2), respectively
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hermite.h>

#include "error.h"
#include "eval.h"

#define pow2(n) (gsl_sf_pow_int(2,n))

/* Evaluates the probabilists' Hermite polynomial of order n at position x using upward recurrence. */
static int 
gsl_sf_hermite_prob_iter_e(const int n, const double x, gsl_sf_result * result)
{
  result->val = 0.;
  result->err = 0.;

  if(n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = 1.;
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = x;
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else if(x == 0.){
    if(GSL_IS_ODD(n)){
      result->val = 0.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
    else{
      if(n < 301){
	/*
	double f;
	int j;
	f = (GSL_IS_ODD(n/2)?-1.:1.);
	for(j=1; j < n; j+=2) {
	  f*=j;
	}
	result->val = f;
	result->err = 0.;
	*/
	if(n < 297){
	  gsl_sf_doublefact_e(n-1, result);
	  (GSL_IS_ODD(n/2)?result->val = -result->val:1.);
	}
	else if (n == 298){
	  result->val = (GSL_IS_ODD(n/2)?-1.:1.)*1.25527562259930633890922678431e304;
	  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	}
	else{
	  result->val = (GSL_IS_ODD(n/2)?-1.:1.)*3.7532741115719259533385880851e306;
	  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	}
      }
      else{
	result->val = (GSL_IS_ODD(n/2)?GSL_NEGINF:GSL_POSINF);
	result->err = GSL_POSINF;
      }
      return GSL_SUCCESS;
    }
  }
/*
  else if(x*x < 4.0*n && n > 100000) {
    // asymptotic formula
    double f = 1.0;
    int j;
    if(GSL_IS_ODD(n)) {
      f=gsl_sf_fact((n-1)/2)*gsl_sf_pow_int(2,n/2)*M_SQRT2/M_SQRTPI;
    }
    else {
      for(j=1; j < n; j+=2) {
	f*=j;
      }
    }
    return f*exp(x*x/4)*cos(x*sqrt(n)-(n%4)*M_PI_2)/sqrt(sqrt(1-x*x/4.0/n));
    // return f*exp(x*x/4)*cos(x*sqrt(n)-n*M_PI_2)/sqrt(sqrt(1-x*x/4.0/n));
  }
*/
  else{
    /* upward recurrence: He_{n+1} = x He_n - n He_{n-1} */

    double p_n0 = 1.0;  /* He_0(x) */
    double p_n1 = x;    /* He_1(x) */
    double p_n = p_n1;

    double e_n0 = GSL_DBL_EPSILON;
    double e_n1 = fabs(x)*GSL_DBL_EPSILON;
    double e_n = e_n1;

    int j=0, c=0;

    for(j=1; j <= n-1; j++){
      if (gsl_isnan(p_n) == 1){
	break;
      }
      p_n  = x*p_n1-j*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;

      e_n  = (fabs(x)*e_n1+j*e_n0);
      e_n0 = e_n1;
      e_n1 = e_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	e_n0 *= 0.5;
	e_n1 *= 0.5;
	e_n = e_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	e_n0 *= 2.0;
	e_n1 *= 2.0;
	e_n = e_n1;
	c--;
      }
    }

    /*
    // check to see that the correct values are computed, even when overflow strikes in the end; works, thus very large results are accessible by determining mantissa and exponent separately
    double lg2 = 0.30102999566398119521467838;
    double ln10 = 2.3025850929940456840179914546843642076011014886;
    printf("res= %g\n", p_n*pow(10.,((lg2*c)-((long)(lg2*c)))) );
    printf("res= %g * 10^(%ld)\n", p_n*pow(10.,((lg2*c)-((long)(lg2*c))))/pow(10.,((long)(log(fabs(p_n*pow(10.,((lg2*c)-((long)(lg2*c))))))/ln10))), ((long)(log(fabs(p_n*pow(10.,((lg2*c)-((long)(lg2*c))))))/ln10))+((long)(lg2*c)) );
    */

    result->val = pow2(c)*p_n;
    result->err = pow2(c)*e_n + fabs(result->val)*GSL_DBL_EPSILON;
    /* result->err = e_n + n*fabs(p_n)*GSL_DBL_EPSILON;
       no idea, where the factor n came from => removed
     */

    if (gsl_isnan(result->val) != 1){
      return GSL_SUCCESS;
    }
    else{
      return GSL_ERANGE;
    }
  }
}

/* Approximatively evaluates the probabilists' Hermite polynomial of order n at position x.
 * An approximation depending on the x-range (see Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society) is used. */
static int 
gsl_sf_hermite_prob_appr_e(const int n, const double x, gsl_sf_result * result)
{
  /* Plancherel-Rotach approximation (note: Szego defines the Airy function differently!) */
  const double aizero1 = -2.3381074104597670384891972524467; /* first zero of the Airy function Ai */
  double z = fabs(x)*M_SQRT1_2;
  double f = 1.;
  int j;
  for(j=1; j <= n; j++) {
    f*=sqrt(j);
  }
  if (z < sqrt(2*n+1.)+aizero1/M_SQRT2/pow(n,1/6.)){
    double phi = acos(z/sqrt(2*n+1.));
    result->val = f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(2./n,0.25)/sqrt(M_SQRTPI*sin(phi))*sin(M_PI*0.75+(0.5*n+0.25)*(sin(2*phi)-2*phi))*exp(0.5*z*z);
    result->err = 2. * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if (z > sqrt(2*n+1.)-aizero1/M_SQRT2/pow(n,1/6.)){
    double phi = acosh(z/sqrt(2*n+1.));
    result->val = f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(n,-0.25)/M_SQRT2/sqrt(M_SQRT2*M_SQRTPI*sinh(phi))*exp((0.5*n+0.25)*(2*phi-sinh(2*phi)))*exp(0.5*z*z);
    result->err = 2. * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else{
    gsl_sf_result Ai;
    gsl_sf_airy_Ai_e((z-sqrt(2*n+1.))*M_SQRT2*pow(n,1/6.),0,&Ai);
    result->val = f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*sqrt(M_SQRTPI*M_SQRT2)*pow(n,-1/12.)*Ai.val*exp(0.5*z*z);
    result->err = f*sqrt(M_SQRTPI*M_SQRT2)*pow(n,-1/12.)*Ai.err*exp(0.5*z*z) + GSL_DBL_EPSILON*fabs(result->val);
    return GSL_SUCCESS;
  }
}

/* Evaluates the probabilists' Hermite polynomial of order n at position x.
 * For small n upward recurrence is employed, while for large n and NaNs from the iteration an approximation depending on the x-range (see Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society) is used. */
int 
gsl_sf_hermite_prob_e(const int n, const double x, gsl_sf_result * result)
{
  if( (x==0. || n<=100000) && (gsl_sf_hermite_prob_iter_e(n,x,result)==GSL_SUCCESS) ){
    return GSL_SUCCESS;
  }
  else{
    return gsl_sf_hermite_prob_appr_e(n,x,result);
  }
}

double gsl_sf_hermite_prob(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_prob_e(n, x, &result));
}

/* Evaluates the m-th derivative of the probabilists' Hermite polynomial of order n at position x.
 * The direct formula He^{(m)}_n = n!/(n-m)!*He_{n-m}(x) (where He_j(x) is the j-th probabilists' Hermite polynomial and He^{(m)}_j(x) its m-th derivative) is employed. */
int
gsl_sf_hermite_prob_der_e(const int m, const int n, const double x, gsl_sf_result * result)
{
  if(n < 0 || m < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n < m) {
    result->val = 0.;
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else{
    double f = gsl_sf_choose(n,m)*gsl_sf_fact(m);
    gsl_sf_result He;
    gsl_sf_hermite_prob_e(n-m,x,&He);
    result->val = He.val*f;
    result->err = He.err*f + GSL_DBL_EPSILON*fabs(result->val);
    return GSL_SUCCESS;
  }
}

double
gsl_sf_hermite_prob_der(const int m, const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_prob_der_e(m, n, x, &result));
}

/* Evaluates the physicists' Hermite polynomial of order n at position x.
 * For small n upward recurrence is employed, while for large n and NaNs from the iteration an approximation depending on the x-range (see Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society) is used. */
int
gsl_sf_hermite_phys_e(const int n, const double x, gsl_sf_result * result)
{
  result->val = 0.;
  result->err = 0.;

  if(n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = 1.;
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = 2.0*x;
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else if(x == 0.){
    if(GSL_IS_ODD(n)){
      result->val = 0.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
    else{
      if(n < 269){
	double f = pow2(n/2);
	gsl_sf_doublefact_e(n-1, result);
	result->val *= f;
	result->err *= f;
	(GSL_IS_ODD(n/2)?result->val = -result->val:1.);
	/*
	double f;
	int j;
	f = (GSL_IS_ODD(n/2)?-1.:1.);
	for(j=1; j < n; j+=2) {
	  f*=2*j;
	}
	result->val = f;
	result->err = 0.;
	*/
      }
      else{
	result->val = (GSL_IS_ODD(n/2)?GSL_NEGINF:GSL_POSINF);
	result->err = GSL_POSINF;
      }
      return GSL_SUCCESS;
    }
  }
  /*
  else if(x*x < 2.0*n && n > 100000) {
    // asymptotic formula
    double f = 1.0;
    int j;
    if(GSL_IS_ODD(n)) {
      f=gsl_sf_fact((n-1)/2)*gsl_sf_pow_int(2,n)/M_SQRTPI;
    }
    else {
      for(j=1; j < n; j+=2) {
	f*=j;
      }
      f*=gsl_sf_pow_int(2,n/2);
    }
    return f*exp(x*x/2)*cos(x*sqrt(2.0*n)-(n%4)*M_PI_2)/sqrt(sqrt(1-x*x/2.0/n));
    // return f*exp(x*x/2)*cos(x*sqrt(2.0*n)-n*M_PI_2)/sqrt(sqrt(1-x*x/2.0/n));
  }
  */
  else if (n <= 100000){
    /* upward recurrence: H_{n+1} = 2x H_n - 2j H_{n-1} */

    double p_n0 = 1.0;    /* H_0(x) */
    double p_n1 = 2.0*x;  /* H_1(x) */
    double p_n = p_n1;

    double e_n0 = GSL_DBL_EPSILON;
    double e_n1 = 2.*fabs(x)*GSL_DBL_EPSILON;
    double e_n = e_n1;

    int j=0, c=0;

    for(j=1; j <= n-1; j++){
      if (gsl_isnan(p_n) == 1){
	break;
      }
      p_n  = 2.0*x*p_n1-2.0*j*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;

      e_n  = 2.*(fabs(x)*e_n1+j*e_n0);
      e_n0 = e_n1;
      e_n1 = e_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	e_n0 *= 0.5;
	e_n1 *= 0.5;
	e_n = e_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	e_n0 *= 2.0;
	e_n1 *= 2.0;
	e_n = e_n1;
	c--;
      }
    }

    result->val = pow2(c)*p_n;
    result->err = pow2(c)*e_n + fabs(result->val)*GSL_DBL_EPSILON;
    /* result->err = e_n + n*fabs(p_n)*GSL_DBL_EPSILON;
       no idea, where the factor n came from => removed */
    if (gsl_isnan(result->val) != 1){
      return GSL_SUCCESS;
    }
  }

  /* the following condition is implied by the logic above */
  {
    /* Plancherel-Rotach approximation (note: Szego defines the Airy function differently!) */
    const double aizero1 = -2.3381074104597670384891972524467; /* first zero of the Airy function Ai */
    double z = fabs(x);
    double f = 1.;
    int j;
    for(j=1; j <= n; j++) {
      f*=sqrt(j);
    }
    if (z < sqrt(2*n+1.)+aizero1/M_SQRT2/pow(n,1/6.)){
      double phi = acos(z/sqrt(2*n+1.));
      result->val = f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*(GSL_IS_ODD(n)?M_SQRT2:1.)*pow2(n/2)*pow(2./n,0.25)/sqrt(M_SQRTPI*sin(phi))*sin(M_PI*0.75+(0.5*n+0.25)*(sin(2*phi)-2*phi))*exp(0.5*z*z);
      result->err = 2. * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if (z > sqrt(2*n+1.)-aizero1/M_SQRT2/pow(n,1/6.)){
      double phi = acosh(z/sqrt(2*n+1.));
      result->val = f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*(GSL_IS_ODD(n)?1.:M_SQRT1_2)*pow2(n/2)*pow(n,-0.25)/sqrt(M_SQRT2*M_SQRTPI*sinh(phi))*exp((0.5*n+0.25)*(2*phi-sinh(2*phi)))*exp(0.5*z*z);
      result->err = 2. * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else{
      gsl_sf_result Ai;
      gsl_sf_airy_Ai_e((z-sqrt(2*n+1.))*M_SQRT2*pow(n,1/6.),0,&Ai);
      result->val = f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*(GSL_IS_ODD(n)?M_SQRT2:1.)*sqrt(M_SQRTPI*M_SQRT2)*pow2(n/2)*pow(n,-1/12.)*Ai.val*exp(0.5*z*z);
      result->err = f*(GSL_IS_ODD(n)?M_SQRT2:1.)*pow2(n/2)*sqrt(M_SQRTPI*M_SQRT2)*pow(n,-1/12.)*Ai.err*exp(0.5*z*z) + GSL_DBL_EPSILON*fabs(result->val);
      return GSL_SUCCESS;
    }
  }
}

double
gsl_sf_hermite_phys(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_phys_e(n, x, &result));
}

/* Evaluates the m-th derivative of the physicists' Hermite polynomial of order n at position x.
 * The direct formula H^{(m)}_n = 2**m*n!/(n-m)!*H_{n-m}(x) (where H_j(x) is the j-th physicists' Hermite polynomial and H^{(m)}_j(x) its m-th derivative) is employed. */
int 
gsl_sf_hermite_phys_der_e(const int m, const int n, const double x, gsl_sf_result * result)
{
  if(n < 0 || m < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n < m) {
    result->val = 0.;
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else{
    double f = gsl_sf_choose(n,m)*gsl_sf_fact(m)*pow2(m);
    gsl_sf_result H;
    gsl_sf_hermite_phys_e(n-m,x,&H);
    result->val = H.val*f;
    result->err = H.err*f + GSL_DBL_EPSILON*fabs(result->val);
    return GSL_SUCCESS;
  }
}

double
gsl_sf_hermite_phys_der(const int m, const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_phys_der_e(m, n, x, &result));
}

/* Evaluates the Hermite function of order n at position x.
 * For large n an approximation depending on the x-range (see Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society) is used, while for small n the direct formula via the probabilists' Hermite polynomial is applied. */
int
gsl_sf_hermite_func_e(const int n, const double x, gsl_sf_result * result)
{
  /*
  if (x*x < 2.0*n && n > 100000){
    // asymptotic formula
    double f = 1.0;
    int j;
    // return f*exp(x*x/4)*cos(x*sqrt(n)-n*M_PI_2)/sqrt(sqrt(1-x*x/4.0/n));
    return cos(x*sqrt(2.0*n)-(n%4)*M_PI_2)/sqrt(sqrt(0.5*n/M_PI*(1-0.5*x*x/n)))/M_PI;
  }
  */
  if (n < 0){
    DOMAIN_ERROR(result);
  }
  else if(n == 0 && x != 0.) {
    result->val = exp(-0.5*x*x)/sqrt(M_SQRTPI);
    result->err = GSL_DBL_EPSILON*fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(n == 1 && x != 0.) {
    result->val = M_SQRT2*x*exp(-0.5*x*x)/sqrt(M_SQRTPI);
    result->err = GSL_DBL_EPSILON*fabs(result->val);
    return GSL_SUCCESS;
  }
  else if (x == 0.){
    if (GSL_IS_ODD(n)){
      result->val = 0.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
    else{
      double f;
      int j;
      f = (GSL_IS_ODD(n/2)?-1.:1.);
      for(j=1; j < n; j+=2) {
	f*=sqrt(j/(j+1.));
      }
      result->val = f/sqrt(M_SQRTPI);
      result->err = GSL_DBL_EPSILON*fabs(result->val);
      return GSL_SUCCESS;
    }
  }
  else if (n <= 100000){
    double f = exp(-0.5*x*x)/sqrt(M_SQRTPI*gsl_sf_fact(n));
    gsl_sf_result He;
    gsl_sf_hermite_prob_iter_e(n,M_SQRT2*x,&He);
    result->val = He.val*f;
    result->err = He.err*f + GSL_DBL_EPSILON*fabs(result->val);
    if (gsl_isnan(result->val) != 1 && f > GSL_DBL_MIN && gsl_finite(He.val) == 1){
      return GSL_SUCCESS;
    }
  }

  /* upward recurrence: Psi_{n+1} = sqrt(2/(n+1))*x Psi_n - sqrt(n/(n+1)) Psi_{n-1} */

  {
    double tw = exp(-x*x*0.5/n); /* "twiddle factor" (in the spirit of FFT) */
    double p_n0 = tw/sqrt(M_SQRTPI); /* Psi_0(x) */
    double p_n1 = p_n0*M_SQRT2*x;    /* Psi_1(x) */
    double p_n = p_n1;
    double e_n0 = p_n0*GSL_DBL_EPSILON;
    double e_n1 = p_n1*GSL_DBL_EPSILON;
    double e_n = e_n1;
  
    int j;

    int c = 0;
    for (j=1;j<n;j++)
      {
        if (gsl_isnan(p_n) == 1){
	  break;
        }
        p_n=tw*(M_SQRT2*x*p_n1-sqrt(j)*p_n0)/sqrt(j+1.);
        p_n0=tw*p_n1;
        p_n1=p_n;

        e_n=tw*(M_SQRT2*fabs(x)*e_n1+sqrt(j)*e_n0)/sqrt(j+1.);
        e_n0=tw*e_n1;
        e_n1=e_n;

        while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
  	p_n0 *= 0.5;
  	p_n1 *= 0.5;
	p_n = p_n1;
	e_n0 *= 0.5;
	e_n1 *= 0.5;
	e_n = e_n1;
	c++;
      }

	while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 = p_n0*2;
	p_n1 = p_n1*2;
	p_n = p_n1;
	e_n0 = e_n0*2;
	e_n1 = e_n1*2;
	e_n = e_n1;
	c--;
      }
    }

  result->val = p_n*pow2(c);
  result->err = n*fabs(result->val)*GSL_DBL_EPSILON;

  if (gsl_isnan(result->val) != 1){
    return GSL_SUCCESS;
  }

  {
    /* Plancherel-Rotach approximation (note: Szego defines the Airy function differently!) */
    const double aizero1 = -2.3381074104597670384891972524467; /* first zero of the Airy function Ai */
    double z = fabs(x);
    if (z < sqrt(2*n+1.)+aizero1/M_SQRT2/pow(n,1/6.)){
      double phi = acos(z/sqrt(2*n+1.));
      result->val = (GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(2./n,0.25)/M_SQRTPI/sqrt(sin(phi))*sin(M_PI*0.75+(0.5*n+0.25)*(sin(2*phi)-2*phi));
      result->err = 2. * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if (z > sqrt(2*n+1.)-aizero1/M_SQRT2/pow(n,1/6.)){
      double phi = acosh(z/sqrt(2*n+1.));
      result->val = (GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(n,-0.25)/
2/M_SQRTPI/sqrt(sinh(phi)/M_SQRT2)*exp((0.5*n+0.25)*(2*phi-sinh(2*phi)));
      result->err = 2. * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else{
      gsl_sf_result Ai;
      gsl_sf_airy_Ai_e((z-sqrt(2*n+1.))*M_SQRT2*pow(n,1/6.),0,&Ai);
      result->val = (GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*sqrt(M_SQRT2)*pow(n,-1/12.)*Ai.val;
      result->err = sqrt(M_SQRT2)*pow(n,-1/12.)*Ai.err +  GSL_DBL_EPSILON*fabs(result->val);
      return GSL_SUCCESS;
    }
  }
  }
}

double
gsl_sf_hermite_func(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_func_e(n, x, &result));
}

/* Evaluates all probabilists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_prob_array(const int nmax, const double x, double * result_array)
{
  if(nmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(nmax == 0) {
    result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else if(nmax == 1) {
    result_array[0] = 1.0;
    result_array[1] = x;
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence: He_{n+1} = x He_n - n He_{n-1} */

    double p_n0 = 1.0;    /* He_0(x) */
    double p_n1 = x;      /* He_1(x) */
    double p_n = p_n1;
    int j=0, c=0;

    result_array[0] = 1.0;
    result_array[1] = x;

    for(j=1; j <= nmax-1; j++){
      p_n  = x*p_n1-j*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	c--;
      }

      result_array[j+1] = pow2(c)*p_n;
    }

    return GSL_SUCCESS;
  }
}


/* Evaluates the m-th derivative of all probabilists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */

int
gsl_sf_hermite_prob_array_der(const int m, const int nmax, const double x, double * result_array)
{
  if(nmax < 0 || m < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(m == 0) {
    gsl_sf_hermite_prob_array(nmax, x, result_array);
    return GSL_SUCCESS;
  }
  else if(nmax < m) {
    int j;
    for(j=0; j <= nmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(nmax == m) {
    int j;
    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }
    result_array[nmax] = gsl_sf_fact(m);
    return GSL_SUCCESS;
  }
  else if(nmax == m+1) {
    int j;
    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }
    result_array[nmax-1] = gsl_sf_fact(m);
    result_array[nmax] = result_array[nmax-1]*(m+1)*x;
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence: He^{(m)}_{n+1} = (n+1)/(n-m+1)*(x He^{(m)}_n - n He^{(m)}_{n-1}) */

    double p_n0 = gsl_sf_fact(m);    /* He^{(m)}_{m}(x) */
    double p_n1 = p_n0*(m+1)*x;      /* He^{(m)}_{m+1}(x) */
    double p_n = p_n1;
    int j=0, c=0;

    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }

    result_array[m] = p_n0;
    result_array[m+1] = p_n1;

    for(j=m+1; j <= nmax-1; j++){
      p_n  = (x*p_n1-j*p_n0)*(j+1)/(j-m+1);
      p_n0 = p_n1;
      p_n1 = p_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	c--;
      }

      result_array[j+1] = pow2(c)*p_n;
    }

    return GSL_SUCCESS;
  }
}

/* Evaluates all derivatives (starting from 0) up to the mmax-th derivative of the probabilists' Hermite polynomial of order n at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */

int
gsl_sf_hermite_prob_der_array(const int mmax, const int n, const double x, double * result_array)
{
  if(n < 0 || mmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 0) {
    int j;
    result_array[0] = 1.0;
    for(j=1; j <= mmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(n == 1 && mmax > 0) {
    int j;
    result_array[0] = x;
    result_array[1] = 1.0;
    for(j=2; j <= mmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if( mmax == 0) {
    result_array[0] = gsl_sf_hermite_prob(n,x);
    return GSL_SUCCESS;
  }
  else if( mmax == 1) {
    result_array[0] = gsl_sf_hermite_prob(n,x);
    result_array[1] = n*gsl_sf_hermite_prob(n-1,x);
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence */

    int k = GSL_MAX_INT(0,n-mmax);
    /* Getting a bit lazy here... */
    double p_n0 = gsl_sf_hermite_prob(k,x);        /* He_k(x) */
    double p_n1 = gsl_sf_hermite_prob(k+1,x);      /* He_{k+1}(x) */
    double p_n  = p_n1;
    int j=0, c=0;

    for(j=n+1; j <= mmax; j++){
      result_array[j] = 0.0;
    }

    result_array[GSL_MIN_INT(n,mmax)] = p_n0;
    result_array[GSL_MIN_INT(n,mmax)-1] = p_n1;

    for(j=GSL_MIN_INT(mmax,n)-1; j > 0; j--){
      k++;
      p_n  = x*p_n1-k*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	c--;
      }

      result_array[j-1] = pow2(c)*p_n;
    }

    p_n = 1.0;
    for(j=1; j <= GSL_MIN_INT(n,mmax); j++){
      p_n  = p_n*(n-j+1);
      result_array[j] = p_n*result_array[j];
    }

    return GSL_SUCCESS;
  }
}

/* Evaluates the series sum_{j=0}^n a_j*He_j(x) with He_j being the j-th probabilists' Hermite polynomial.
 * For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to probabilists' Hermite polynomials is used. */

int
gsl_sf_hermite_prob_series_e(const int n, const double x, const double * a, gsl_sf_result * result)
{
  if(n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = a[0];
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = a[0]+a[1]*x;
    result->err = 2.*GSL_DBL_EPSILON * (fabs(a[0]) + fabs(a[1]*x)) ;
    return GSL_SUCCESS;
  }
  else {
    /* downward recurrence: b_n = a_n + x b_{n+1} - (n+1) b_{n+2} */

    double b0 = 0.;
    double b1 = 0.;
    double btmp = 0.;

    double e0 = 0.;
    double e1 = 0.;
    double etmp = e1;

    int j;

    for(j=n; j >= 0; j--){
      btmp = b0;
      b0  = a[j]+x*b0-(j+1)*b1;
      b1 = btmp;

      etmp = e0;
      e0  = (GSL_DBL_EPSILON*fabs(a[j])+fabs(x)*e0+(j+1)*e1);
      e1 = etmp;
    }

    result->val = b0;
    result->err = e0 + fabs(b0)*GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
}

double
gsl_sf_hermite_prob_series(const int n, const double x, const double * a)
{
  EVAL_RESULT(gsl_sf_hermite_prob_series_e(n, x, a, &result));
}

/* Evaluates all physicists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_phys_array(const int nmax, const double x, double * result_array)
{
  if(nmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(nmax == 0) {
    result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else if(nmax == 1) {
    result_array[0] = 1.0;
    result_array[1] = 2.0*x;
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence: H_{n+1} = 2x H_n - 2n H_{n-1} */

    double p_n0 = 1.0;      /* H_0(x) */
    double p_n1 = 2.0*x;    /* H_1(x) */
    double p_n = p_n1;
    int j=0, c=0;

    result_array[0] = 1.0;
    result_array[1] = 2.0*x;

    for(j=1; j <= nmax-1; j++){
      p_n  = 2.0*x*p_n1-2.0*j*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	c--;
      }

      result_array[j+1] = pow2(c)*p_n;
    }

    return GSL_SUCCESS;
  }
}


/* Evaluates the m-th derivative of all physicists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_phys_array_der(const int m, const int nmax, const double x, double * result_array)
{
  if(nmax < 0 || m < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(m == 0) {
    gsl_sf_hermite_phys_array(nmax, x, result_array);
    return GSL_SUCCESS;
  }
  else if(nmax < m) {
    int j;
    for(j=0; j <= nmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(nmax == m) {
    int j;
    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }
    result_array[nmax] = pow2(m)*gsl_sf_fact(m);
    return GSL_SUCCESS;
  }
  else if(nmax == m+1) {
    int j;
    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }
    result_array[nmax-1] = pow2(m)*gsl_sf_fact(m);
    result_array[nmax] = result_array[nmax-1]*2*(m+1)*x;
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence: H^{(m)}_{n+1} = 2(n+1)/(n-m+1)*(x H^{(m)}_n - n H^{(m)}_{n-1}) */

    double p_n0 = pow2(m)*gsl_sf_fact(m);    /* H^{(m)}_{m}(x) */
    double p_n1 = p_n0*2*(m+1)*x;    /* H^{(m)}_{m+1}(x) */
    double p_n = p_n1;
    int j=0, c=0;

    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }

    result_array[m] = p_n0;
    result_array[m+1] = p_n1;

    for(j=m+1; j <= nmax-1; j++){
      p_n  = (x*p_n1-j*p_n0)*2*(j+1)/(j-m+1);
      p_n0 = p_n1;
      p_n1 = p_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	c--;
      }

      result_array[j+1] = pow2(c)*p_n;
    }

    return GSL_SUCCESS;
  }
}


/* Evaluates all derivatives (starting from 0) up to the mmax-th derivative of the physicists' Hermite polynomial of order n at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_phys_der_array(const int mmax, const int n, const double x, double * result_array)
{
  if(n < 0 || mmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 0) {
    int j;
    result_array[0] = 1.0;
    for(j=1; j <= mmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(n == 1 && mmax > 0) {
    int j;
    result_array[0] = 2*x;
    result_array[1] = 1.0;
    for(j=2; j <= mmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if( mmax == 0) {
    result_array[0] = gsl_sf_hermite_phys(n,x);
    return GSL_SUCCESS;
  }
  else if( mmax == 1) {
    result_array[0] = gsl_sf_hermite_phys(n,x);
    result_array[1] = 2*n*gsl_sf_hermite_phys(n-1,x);
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence */

    int k = GSL_MAX_INT(0,n-mmax);
    /* Getting a bit lazy here... */
    double p_n0 = gsl_sf_hermite_phys(k,x);        /* H_k(x) */
    double p_n1 = gsl_sf_hermite_phys(k+1,x);      /* H_{k+1}(x) */
    double p_n  = p_n1;
    int j=0, c=0;

    for(j=n+1; j <= mmax; j++){
      result_array[j] = 0.0;
    }

    result_array[GSL_MIN_INT(n,mmax)] = p_n0;
    result_array[GSL_MIN_INT(n,mmax)-1] = p_n1;

    for(j=GSL_MIN_INT(mmax,n)-1; j > 0; j--){
      k++;
      p_n  = 2*x*p_n1-2*k*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	c--;
      }

      result_array[j-1] = pow2(c)*p_n;
    }

    p_n = 1.0;
    for(j=1; j <= GSL_MIN_INT(n,mmax); j++){
      p_n  = p_n*(n-j+1)*2;
      result_array[j] = p_n*result_array[j];
    }

    return GSL_SUCCESS;
  }
}


/* Evaluates the series sum_{j=0}^n a_j*H_j(x) with H_j being the j-th physicists' Hermite polynomial.
 * For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to physicists' Hermite polynomials is used. */
int
gsl_sf_hermite_phys_series_e(const int n, const double x, const double * a, gsl_sf_result * result)
{
  if(n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = a[0];
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = a[0]+a[1]*2.*x;
    result->err = 2.*GSL_DBL_EPSILON * (fabs(a[0]) + fabs(a[1]*2.*x)) ;
    return GSL_SUCCESS;
  }
  else {
    /* downward recurrence: b_n = a_n + 2x b_{n+1} - 2(n+1) b_{n+2} */

    double b0 = 0.;
    double b1 = 0.;
    double btmp = 0.;

    double e0 = 0.;
    double e1 = 0.;
    double etmp = e1;

    int j;

    for(j=n; j >= 0; j--){
      btmp = b0;
      b0  = a[j]+2.*x*b0-2.*(j+1)*b1;
      b1 = btmp;

      etmp = e0;
      e0  = (GSL_DBL_EPSILON*fabs(a[j])+fabs(2.*x)*e0+2.*(j+1)*e1);
      e1 = etmp;
    }

    result->val = b0;
    result->err = e0 + fabs(b0)*GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
}

double
gsl_sf_hermite_phys_series(const int n, const double x, const double * a)
{
  EVAL_RESULT(gsl_sf_hermite_phys_series_e(n, x, a, &result));
}


/* Evaluates all Hermite functions up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_func_array(const int nmax, const double x, double * result_array)
{
  if(nmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(nmax == 0) {
    result_array[0] = exp(-0.5*x*x)/sqrt(M_SQRTPI);
    return GSL_SUCCESS;
  }
  else if(nmax == 1) {
    result_array[0] = exp(-0.5*x*x)/sqrt(M_SQRTPI);
    result_array[1] = result_array[0]*M_SQRT2*x;
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence: Psi_{n+1} = sqrt(2/(n+1))*x Psi_n - sqrt(n/(n+1)) Psi_{n-1} */

    double p_n0 = exp(-0.5*x*x)/sqrt(M_SQRTPI);   /* Psi_0(x) */
    double p_n1 = p_n0*M_SQRT2*x;                 /* Psi_1(x) */
    double p_n = p_n1;
    int j=0, c=0;

    result_array[0] = p_n0;
    result_array[1] = p_n1;

  for (j=1;j<=nmax-1;j++)
    {
      p_n=(M_SQRT2*x*p_n1-sqrt(j)*p_n0)/sqrt(j+1.);
      p_n0=p_n1;
      p_n1=p_n;

      while(( GSL_MIN(fabs(p_n0),fabs(p_n1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) > GSL_SQRT_DBL_MAX )){
	p_n0 *= 0.5;
	p_n1 *= 0.5;
	p_n = p_n1;
	c++;
      }

      while(( ( ( fabs(p_n0) < GSL_SQRT_DBL_MIN ) && ( p_n0 != 0) ) || ( ( fabs(p_n1) < GSL_SQRT_DBL_MIN ) && ( p_n1 != 0) ) ) && ( GSL_MAX(fabs(p_n0),fabs(p_n1)) < 0.5*GSL_SQRT_DBL_MAX )){
	p_n0 *= 2.0;
	p_n1 *= 2.0;
	p_n = p_n1;
	c--;
      }

      result_array[j+1] = pow2(c)*p_n;
    }

    return GSL_SUCCESS;
  }
}

/* Evaluates the series sum_{j=0}^n a_j*Psi_j(x) with Psi_j being the j-th Hermite function.
 * For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to Hermite functions is used. */

int
gsl_sf_hermite_func_series_e(const int n, const double x, const double * a, gsl_sf_result * result)
{
  if(n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = a[0]*exp(-0.5*x*x)/sqrt(M_SQRTPI);
    result->err = GSL_DBL_EPSILON*fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = (a[0]+a[1]*M_SQRT2*x)*exp(-0.5*x*x)/sqrt(M_SQRTPI);
    result->err = 2.*GSL_DBL_EPSILON*(fabs(a[0])+fabs(a[1]*M_SQRT2*x))*exp(-0.5*x*x)/sqrt(M_SQRTPI);
    return GSL_SUCCESS;
  }
  else {
    /* downward recurrence: b_n = a_n + sqrt(2/(n+1))*x b_{n+1} - sqrt((n+1)/(n+2)) b_{n+2} */

    double b0 = 0.;
    double b1 = 0.;
    double btmp = 0.;

    double e0 = 0.;
    double e1 = 0.;
    double etmp = e1;

    int j;

    for(j=n; j >= 0; j--){
      btmp = b0;
      b0  = a[j]+sqrt(2./(j+1))*x*b0-sqrt((j+1.)/(j+2.))*b1;
      b1 = btmp;

      etmp = e0;
      e0  = (GSL_DBL_EPSILON*fabs(a[j])+sqrt(2./(j+1))*fabs(x)*e0+sqrt((j+1.)/(j+2.))*e1);
      e1 = etmp;
    }

    result->val = b0*exp(-0.5*x*x)/sqrt(M_SQRTPI);
    result->err = e0 + fabs(result->val)*GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
}

double
gsl_sf_hermite_func_series(const int n, const double x, const double * a)
{
  EVAL_RESULT(gsl_sf_hermite_func_series_e(n, x, a, &result));
}


/* Evaluates the m-th derivative of the Hermite function of order n at position x.
 * A summation including upward recurrences is used. */
int
gsl_sf_hermite_func_der_e(const int m, const int n, const double x, gsl_sf_result * result)
{
  /* FIXME: asymptotic formula! */
  if(m < 0 || n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(m == 0){
    return gsl_sf_hermite_func_e(n,x,result);
  }
  else{
    int j=0, c=0;
    double r,er,b;
    double h0 = 1.;
    double h1 = x;
    double eh0 = GSL_DBL_EPSILON;
    double eh1 = GSL_DBL_EPSILON;
    double p0 = 1.;
    double p1 = M_SQRT2*x;
    double ep0 = GSL_DBL_EPSILON;
    double ep1 = M_SQRT2*GSL_DBL_EPSILON;
    double f = 1.;
    for (j=GSL_MAX_INT(1,n-m+1);j<=n;j++)
      {
	f *= sqrt(2.*j);
      }
    if (m>n)
      {
	f = (GSL_IS_ODD(m-n)?-f:f);
	for (j=0;j<GSL_MIN_INT(n,m-n);j++)
	  {
	    f *= (m-j)/(j+1.);
	  }
      }
    c = 0;
    for (j=1;j<=m-n;j++)
      {
	b = x*h1-j*h0;
	h0 = h1;
	h1 = b;
	
	b  = (fabs(x)*eh1+j*eh0);
	eh0 = eh1;
	eh1 = b;

	while(( GSL_MIN(fabs(h0),fabs(h1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(h0),fabs(h1)) > GSL_SQRT_DBL_MAX )){
	  h0 *= 0.5;
	  h1 *= 0.5;
	  eh0 *= 0.5;
	  eh1 *= 0.5;
	  c++;
	}

	while(( ( (fabs(h0) < GSL_SQRT_DBL_MIN) && (h0 != 0) ) || ( (fabs(h1) < GSL_SQRT_DBL_MIN) && (h1 != 0) ) ) && ( GSL_MAX(fabs(h0),fabs(h1)) < 0.5*GSL_SQRT_DBL_MAX )){
	  h0 *= 2.0;
	  h1 *= 2.0;
	  eh0 *= 2.0;
	  eh1 *= 2.0;
	  c--;
	}

      }
    h0 *= pow2(c);
    h1 *= pow2(c);
    eh0 *= pow2(c);
    eh1 *= pow2(c);

    b = 0.;
    c = 0;
    for (j=1;j<=n-m;j++)
      {
	b = (M_SQRT2*x*p1-sqrt(j)*p0)/sqrt(j+1.);
	p0 = p1;
	p1 = b;
	
	b  = (M_SQRT2*fabs(x)*ep1+sqrt(j)*ep0)/sqrt(j+1.);
	ep0 = ep1;
	ep1 = b;

	while(( GSL_MIN(fabs(p0),fabs(p1)) > 2.0*GSL_SQRT_DBL_MIN ) && ( GSL_MAX(fabs(p0),fabs(p1)) > GSL_SQRT_DBL_MAX )){
	  p0 *= 0.5;
	  p1 *= 0.5;
	  ep0 *= 0.5;
	  ep1 *= 0.5;
	  c++;
	}

	while(( ( (fabs(p0) < GSL_SQRT_DBL_MIN) && (p0 != 0) ) || ( (fabs(p1) < GSL_SQRT_DBL_MIN) && (p1 != 0) ) ) && ( GSL_MAX(fabs(p0),fabs(p1)) < 0.5*GSL_SQRT_DBL_MAX )){
	  p0 = p0*2;
	  p1 = p1*2;
	  ep0 = ep0*2;
	  ep1 = ep1*2;
	  c--;
	}

      }
    p0 *= pow2(c);
    p1 *= pow2(c);
    ep0 *= pow2(c);
    ep1 *= pow2(c);

    c = 0;
    b = 0.;
    r = 0.;
    er = 0.;
    for (j=GSL_MAX_INT(0,m-n);j<=m;j++)
      {
	r += f*h0*p0;
	er += eh0*fabs(f*p0)+ep0*fabs(f*h0)+GSL_DBL_EPSILON*fabs(f*h0*p0);
	
	b = x*h1-(j+1.)*h0;
	h0 = h1;
	h1 = b;
	
	b  = 0.5*(fabs(x)*eh1+(j+1.)*eh0);
	eh0 = eh1;
	eh1 = b;
	
	b = (M_SQRT2*x*p1-sqrt(n-m+j+1.)*p0)/sqrt(n-m+j+2.);
	p0 = p1;
	p1 = b;
	
	b  = 0.5*(M_SQRT2*fabs(x)*ep1+sqrt(n-m+j+1.)*ep0)/sqrt(n-m+j+2.);
	ep0 = ep1;
	ep1 = b;
	
	f *= -(m-j)/(j+1.)/sqrt(n-m+j+1.)*M_SQRT1_2;

	while(( (fabs(h0) > 2.0*GSL_SQRT_DBL_MIN) && (fabs(h1) > 2.0*GSL_SQRT_DBL_MIN) && (fabs(p0) > 2.0*GSL_SQRT_DBL_MIN) && (fabs(p1) > 2.0*GSL_SQRT_DBL_MIN) && (fabs(r) > 4.0*GSL_SQRT_DBL_MIN) ) && ( (fabs(h0) > GSL_SQRT_DBL_MAX) || (fabs(h1) > GSL_SQRT_DBL_MAX) || (fabs(p0) > GSL_SQRT_DBL_MAX) || (fabs(p1) > GSL_SQRT_DBL_MAX) || (fabs(r) > GSL_SQRT_DBL_MAX) )){
	  h0 *= 0.5;
	  h1 *= 0.5;
	  eh0 *= 0.5;
	  eh1 *= 0.5;
	  p0 *= 0.5;
	  p1 *= 0.5;
	  ep0 *= 0.5;
	  ep1 *= 0.5;
	  r *= 0.25;
	  er *= 0.25;
	  c++;
	}

	while(( ( (fabs(h0) < GSL_SQRT_DBL_MIN) && (h0 != 0) ) || ( (fabs(h1) < GSL_SQRT_DBL_MIN) && (h1 != 0) ) || ( (fabs(p0) < GSL_SQRT_DBL_MIN) && (p0 != 0) ) || ( (fabs(p1) < GSL_SQRT_DBL_MIN) && (p1 != 0) ) || ( (fabs(r) < GSL_SQRT_DBL_MIN) && (r != 0) ) ) && ( (fabs(h0) < 0.5*GSL_SQRT_DBL_MAX) && (fabs(h1) < 0.5*GSL_SQRT_DBL_MAX) && (fabs(p0) < 0.5*GSL_SQRT_DBL_MAX) && (fabs(p1) < 0.5*GSL_SQRT_DBL_MAX) && (fabs(r) < 0.25*GSL_SQRT_DBL_MAX) )){
	  p0 *= 2.0;
	  p1 *= 2.0;
	  ep0 *= 2.0;
	  ep1 *= 2.0;
	  h0 *= 2.0;
	  h1 *= 2.0;
	  eh0 *= 2.0;
	  eh1 *= 2.0;
	  r *= 4.0;
	  er *= 4.0;
	  c--;
	}

      }

    r *= pow2(2*c);
    er *= pow2(2*c);
    result->val = r*exp(-0.5*x*x)/sqrt(M_SQRTPI);
    result->err = er*fabs(exp(-0.5*x*x)/sqrt(M_SQRTPI)) + GSL_DBL_EPSILON*fabs(result->val);
    return GSL_SUCCESS;
  }
}

double
gsl_sf_hermite_func_der(const int m, const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_func_der_e(m, n, x, &result));
}

static double
H_zero_init(const int n, const int k)
{
  double p = 1., x = 1., y = 1.;
  if (k == 1 && n > 50) {
    x = (GSL_IS_ODD(n)?1./sqrt((n-1)/6.):1./sqrt(0.5*n));
  }
  else {
    p = -0.7937005259840997373758528196*gsl_sf_airy_zero_Ai(n/2-k+1);
    x = sqrt(2*n+1.);
    y = pow(2*n+1.,1/6.);
    x = x - p/y - 0.1*p*p/(x*y*y) + (9/280. - p*p*p*11/350.)/(x*x*x) + (p*277/12600. - gsl_sf_pow_int(p,4)*823/63000.)/gsl_sf_pow_int(x,4)/y;
  }
  p = acos(x/sqrt(2*n+1.));
  y = M_PI*(-2*(n/2-k)-1.5)/(n+0.5);
  if(gsl_fcmp(y,sin(2.*p)-2*p,GSL_SQRT_DBL_EPSILON)==0) return x; /* initial approx sufficiently accurate */
  if (y > -GSL_DBL_EPSILON) return sqrt(2*n+1.);
  if (p < GSL_DBL_EPSILON) p = GSL_DBL_EPSILON;
  if (p > M_PI_2) p = M_PI_2;
  if (sin(2.*p)-2*p > y){
    x = GSL_MAX((sin(2.*p)-2*p-y)/4.,GSL_SQRT_DBL_EPSILON);
    do{
      x *= 2.;
      p += x;
    } while (sin(2.*p)-2*p > y);
  }
  do {
    x = p;
    p -= (sin(2.*p)-2.*p-y)/(2.*cos(2.*p)-2.);
    if (p<0.||p>M_PI_2) p = M_PI_2;
  } while (gsl_fcmp(x,p,100*GSL_DBL_EPSILON)!=0);
  return sqrt(2*n+1.)*cos(p);
}


/* lookup table for the positive zeros of the probabilists' Hermite polynomials of order 3 through 20 */
static double He_zero_tab[99] = {
  1.73205080756887729352744634151,
  0.741963784302725857648513596726,
  2.33441421833897723931751226721,
  1.35562617997426586583052129087,
  2.85697001387280565416230426401,
  0.616706590192594152193686099399,
  1.88917587775371067550566789858,
  3.32425743355211895236183546247,
  1.154405394739968127239597758838,
  2.36675941073454128861885646856,
  3.75043971772574225630392202571,
  0.539079811351375108072461918694,
  1.63651904243510799922544657297,
  2.80248586128754169911301080618,
  4.14454718612589433206019783917,
  1.023255663789132524828148225810,
  2.07684797867783010652215614374,
  3.20542900285646994336567590292,
  4.51274586339978266756667884317,
  0.484935707515497653046233483105,
  1.46598909439115818325066466416,
  2.48432584163895458087625118368,
  3.58182348355192692277623675546,
  4.85946282833231215015516494660,
  0.928868997381063940144111999584,
  1.87603502015484584534137013967,
  2.86512316064364499771968407254,
  3.93616660712997692868589612142,
  5.18800122437487094818666404539,
  0.444403001944138945299732445510,
  1.34037519715161672153112945211,
  2.25946445100079912386492979448,
  3.22370982877009747166319001956,
  4.27182584793228172295999293076,
  5.50090170446774760081221630899,
  0.856679493519450033897376121795,
  1.72541837958823916151095838741,
  2.62068997343221478063807762201,
  3.56344438028163409162493844661,
  4.59139844893652062705231872720,
  5.80016725238650030586450565322,
  0.412590457954601838167454145167,
  1.24268895548546417895063983219,
  2.08834474570194417097139675101,
  2.96303657983866750254927123447,
  3.88692457505976938384755016476,
  4.89693639734556468372449782879,
  6.08740954690129132226890147034,
  0.799129068324547999424888414207,
  1.60671006902872973652322479373,
  2.43243682700975804116311571682,
  3.28908242439876638890856229770,
  4.19620771126901565957404160583,
  5.19009359130478119946445431715,
  6.36394788882983831771116094427,
  0.386760604500557347721047189801,
  1.16382910055496477419336819907,
  1.95198034571633346449212362880,
  2.76024504763070161684598142269,
  3.60087362417154828824902745506,
  4.49295530252001124266582263095,
  5.47222570594934308841242925805,
  6.63087819839312848022981922233,
  0.751842600703896170737870774614,
  1.50988330779674075905491513417,
  2.28101944025298889535537879396,
  3.07379717532819355851658337833,
  3.90006571719800990903311840097,
  4.77853158962998382710540812497,
  5.74446007865940618125547815768,
  6.88912243989533223256205432938,
  0.365245755507697595916901619097,
  1.09839551809150122773848360538,
  1.83977992150864548966395498992,
  2.59583368891124032910545091458,
  3.37473653577809099529779309480,
  4.18802023162940370448450911428,
  5.05407268544273984538327527397,
  6.00774591135959752029303858752,
  7.13946484914647887560975631213,
  0.712085044042379940413609979021,
  1.42887667607837287134157901452,
  2.15550276131693514033871248449,
  2.89805127651575312007902775275,
  3.66441654745063847665304033851,
  4.46587262683103133615452574019,
  5.32053637733603803162823765939,
  6.26289115651325170419416064557,
  7.38257902403043186766326977122,
  0.346964157081355927973322447164,
  1.04294534880275103146136681143,
  1.74524732081412671493067861704,
  2.45866361117236775131735057433,
  3.18901481655338941485371744116,
  3.94396735065731626033176813604,
  4.73458133404605534390170946748,
  5.57873880589320115268040332802,
  6.51059015701365448636289263918,
  7.61904854167975829138128156060
};

/* Computes the s-th zero the probabilists' Hermite polynomial of order n.
A Newton iteration using a continued fraction representation adapted from [E.T. Whittaker (1914), On the continued fractions which represent the functions of Hermite and other functions defined by differential equations, Proceedings of the Edinburgh Mathematical Society, 32, 65-74] is performed with the initial approximation from [Arpad Elbert and Martin E. Muldoon, Approximations for zeros of Hermite functions, pp. 117-126 in D. Dominici and R. S. Maier, eds, "Special Functions and Orthogonal Polynomials", Contemporary Mathematics, vol 471 (2008)] refined via the bisection method. */
int
gsl_sf_hermite_prob_zero_e(const int n, const int s, gsl_sf_result * result)
{
  if(n <= 0 || s < 0 || s > n/2) {
    DOMAIN_ERROR(result);
  }
  else if(s == 0) {
    if (GSL_IS_ODD(n) == 1) {
      result->val = 0.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
    else {
      DOMAIN_ERROR(result);
    }
  }
  else if(n == 2) {
    result->val = 1.;
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else if(n < 21) {
    result->val = He_zero_tab[(GSL_IS_ODD(n)?n/2:0)+((n/2)*(n/2-1))+s-2];
    result->err = GSL_DBL_EPSILON*(result->val);
    return GSL_SUCCESS;
  }
  else {
    double d = 1., x = 1., x0 = 1.;
    int j;
    x = H_zero_init(n,s) * M_SQRT2;
    do {
      x0 = x;
      d = 0.;
      for (j=1; j<n; j++) d = j/(x-d);
      x -= (x-d)/n;
    /* gsl_fcmp can be used since the smallest zero approaches 1/sqrt(n) or 1/sqrt((n-1)/3.) for large n and thus all zeros are non-zero (except for the trivial case handled above) */
    } while (gsl_fcmp(x,x0,10*GSL_DBL_EPSILON)!=0);
    result->val = x;
    result->err = 2*GSL_DBL_EPSILON*x + fabs(x-x0);
    return GSL_SUCCESS;
  }
}

double
gsl_sf_hermite_prob_zero(const int n, const int s)
{
  EVAL_RESULT(gsl_sf_hermite_prob_zero_e(n, s, &result));
}

/* lookup table for the positive zeros of the physicists' Hermite polynomials of order 3 through 20 */
static double H_zero_tab[99] = {
  1.22474487139158904909864203735,
  0.524647623275290317884060253835,
  1.65068012388578455588334111112,
  0.958572464613818507112770593893,
  2.02018287045608563292872408814,
  0.436077411927616508679215948251,
  1.335849074013696949714895282970,
  2.35060497367449222283392198706,
  0.816287882858964663038710959027,
  1.67355162876747144503180139830,
  2.65196135683523349244708200652,
  0.381186990207322116854718885584,
  1.157193712446780194720765779063,
  1.98165675669584292585463063977,
  2.93063742025724401922350270524,
  0.723551018752837573322639864579,
  1.46855328921666793166701573925,
  2.26658058453184311180209693284,
  3.19099320178152760723004779538,
  0.342901327223704608789165025557,
  1.03661082978951365417749191676,
  1.75668364929988177345140122011,
  2.53273167423278979640896079775,
  3.43615911883773760332672549432,
  0.656809566882099765024611575383,
  1.32655708449493285594973473558,
  2.02594801582575533516591283121,
  2.78329009978165177083671870152,
  3.66847084655958251845837146485,
  0.314240376254359111276611634095,
  0.947788391240163743704578131060,
  1.59768263515260479670966277090,
  2.27950708050105990018772856942,
  3.02063702512088977171067937518,
  3.88972489786978191927164274724,
  0.605763879171060113080537108602,
  1.22005503659074842622205526637,
  1.85310765160151214200350644316,
  2.51973568567823788343040913628,
  3.24660897837240998812205115236,
  4.10133759617863964117891508007,
  0.291745510672562078446113075799,
  0.878713787329399416114679311861,
  1.47668273114114087058350654421,
  2.09518325850771681573497272630,
  2.74847072498540256862499852415,
  3.46265693360227055020891736115,
  4.30444857047363181262129810037,
  0.565069583255575748526020337198,
  1.13611558521092066631913490556,
  1.71999257518648893241583152515,
  2.32573248617385774545404479449,
  2.96716692790560324848896036355,
  3.66995037340445253472922383312,
  4.49999070730939155366438053053,
  0.273481046138152452158280401965,
  0.822951449144655892582454496734,
  1.38025853919888079637208966969,
  1.95178799091625397743465541496,
  2.54620215784748136215932870545,
  3.17699916197995602681399455926,
  3.86944790486012269871942409801,
  4.68873893930581836468849864875,
  0.531633001342654731349086553718,
  1.06764872574345055363045773799,
  1.61292431422123133311288254454,
  2.17350282666662081927537907149,
  2.75776291570388873092640349574,
  3.37893209114149408338327069289,
  4.06194667587547430689245559698,
  4.87134519367440308834927655662,
  0.258267750519096759258116098711,
  0.776682919267411661316659462284,
  1.30092085838961736566626555439,
  1.83553160426162889225383944409,
  2.38629908916668600026459301424,
  2.96137750553160684477863254906,
  3.57376906848626607950067599377,
  4.24811787356812646302342016090,
  5.04836400887446676837203757885,
  0.503520163423888209373811765050,
  1.01036838713431135136859873726,
  1.52417061939353303183354859367,
  2.04923170985061937575050838669,
  2.59113378979454256492128084112,
  3.15784881834760228184318034120,
  3.76218735196402009751489394104,
  4.42853280660377943723498532226,
  5.22027169053748216460967142500,
  0.245340708300901249903836530634,
  0.737473728545394358705605144252,
  1.23407621539532300788581834696,
  1.73853771211658620678086566214,
  2.25497400208927552308233334473,
  2.78880605842813048052503375640,
  3.34785456738321632691492452300,
  3.94476404011562521037562880052,
  4.60368244955074427307767524898,
  5.38748089001123286201690041068
};

/* Computes the s-th zero the physicists' Hermite polynomial of order n, thus also the s-th zero of the Hermite function of order n.
A Newton iteration using a continued fraction representation adapted from [E.T. Whittaker (1914), On the continued fractions which represent the functions of Hermite and other functions defined by differential equations, Proceedings of the Edinburgh Mathematical Society, 32, 65-74] is performed with the initial approximation from [Arpad Elbert and Martin E. Muldoon, Approximations for zeros of Hermite functions, pp. 117-126 in D. Dominici and R. S. Maier, eds, "Special Functions and Orthogonal Polynomials", Contemporary Mathematics, vol 471 (2008)] refined via the bisection method. */
int
gsl_sf_hermite_phys_zero_e(const int n, const int s, gsl_sf_result * result)
{
  if(n <= 0 || s < 0 || s > n/2) {
    DOMAIN_ERROR(result);
  }
  else if(s == 0) {
    if (GSL_IS_ODD(n) == 1) {
      result->val = 0.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
    else {
      DOMAIN_ERROR(result);
    }
  }
  else if(n == 2) {
    result->val = M_SQRT1_2;
    result->err = 0.;
    return GSL_SUCCESS;
  }
  else if(n < 21) {
    result->val = H_zero_tab[(GSL_IS_ODD(n)?n/2:0)+((n/2)*(n/2-1))+s-2];
    result->err = GSL_DBL_EPSILON*(result->val);
    return GSL_SUCCESS;
  }
  else {
    double d = 1., x = 1., x0 = 1.;
    int j;
    x = H_zero_init(n,s);
    do {
      x0 = x;
      d = 0.;
      for (j=1; j<n; j++) d = 2*j/(2.*x-d);
      x -= (2*x-d)*0.5/n;
    /* gsl_fcmp can be used since the smallest zero approaches 1/sqrt(n) or 1/sqrt((n-1)/3.) for large n and thus all zeros are non-zero (except for the trivial case handled above) */
    } while (gsl_fcmp(x,x0,10*GSL_DBL_EPSILON)!=0);
    result->val = x;
    result->err = 2*GSL_DBL_EPSILON*x + fabs(x-x0);
    return GSL_SUCCESS;
  }
}

double
gsl_sf_hermite_phys_zero(const int n, const int s)
{
  EVAL_RESULT(gsl_sf_hermite_phys_zero_e(n, s, &result));
}

int
gsl_sf_hermite_func_zero_e(const int n, const int s, gsl_sf_result * result)
{
  return gsl_sf_hermite_phys_zero_e(n, s, result);
}

double
gsl_sf_hermite_func_zero(const int n, const int s)
{
  EVAL_RESULT(gsl_sf_hermite_func_zero_e(n, s, &result));
}
