/* specfunc/hermite.c
 * 
 * Copyright (C) 2011, 2012, 2013, 2014, 2019 Konrad Griessinger (konradg(at)gmx.net)
 * Copyright (C) 2019 Patrick Alken
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
#define RND(x)  ((double) ((x >= 0) ? (int) (x + 0.5) : (int) (x - 0.5)))

/* evaluates the probabilists' Hermite polynomial of order n at position x */
int 
gsl_sf_hermite_prob_e(const int n, const double x, gsl_sf_result * result)
{
  if (n < 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (n == 0)
    {
      result->val = 1.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
  else if (n == 1)
    {
      result->val = x;
      result->err = 0.;
      return GSL_SUCCESS;
    }
  else if (x == 0.0)
    {
      if (GSL_IS_ODD(n))
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }
      else
        {
          /* for n even, He_n(0) = (-1)^{n/2} (n - 1)!! */

          int status = GSL_SUCCESS;

          if (n - 1 > GSL_SF_DOUBLEFACT_NMAX)  /* test if (n - 1)!! will overflow */
            {
              status = GSL_EOVRFLW;
              result->val = GSL_IS_ODD(n / 2) ? GSL_NEGINF : GSL_POSINF;
              result->err = GSL_POSINF;
            }
          else
            {
              gsl_sf_doublefact_e(n - 1, result);
              if (GSL_IS_ODD(n / 2))
                result->val = -result->val;
            }

          return status;
        }
    }
  else
    {
      /* upward recurrence: He_{n+1} = x He_n - n He_{n-1} */

      int status = GSL_SUCCESS;
      const double abs_x = fabs(x);
      const double thresh1 = abs_x > 1.0 ? 0.9 * GSL_DBL_MAX / abs_x : GSL_DBL_MAX;
      const double thresh2 = 0.9 * GSL_DBL_MAX;

      double p_n0 = 1.0;    /* He_0(x) */
      double p_n1 = x;      /* He_1(x) */
      double p_n = p_n1;

      double e_n0 = GSL_DBL_EPSILON;
      double e_n1 = fabs(x)*GSL_DBL_EPSILON;
      double e_n = e_n1;

      int j;

      for (j = 1; j < n; j++)
        {
          if (fabs(p_n1) > thresh1 ||    /* test if x*p_n1 will overflow */
              fabs(p_n0) > thresh2 / j)  /* test if j*p_n0 will overflow */
            {
              status = GSL_EOVRFLW;
              break;
            }

          p_n  = x*p_n1 - j*p_n0;
          p_n0 = p_n1;
          p_n1 = p_n;

          e_n  = fabs(x)*e_n1+j*e_n0;
          e_n0 = e_n1;
          e_n1 = e_n;
        }

      result->val = p_n;
      result->err = e_n + fabs(result->val)*GSL_DBL_EPSILON;

      return status;
    }
}

double gsl_sf_hermite_prob(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_prob_e(n, x, &result));
}

/* Evaluates the m-th derivative of the probabilists' Hermite polynomial of order n at position x.
 * The direct formula He^{(m)}_n = n!/(n-m)!*He_{n-m}(x) (where He_j(x) is the j-th probabilists' Hermite polynomial and He^{(m)}_j(x) its m-th derivative) is employed. */
int
gsl_sf_hermite_prob_deriv_e(const int m, const int n, const double x, gsl_sf_result * result)
{
  if (n < 0 || m < 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (n < m)
    {
      result->val = 0.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      double f = gsl_sf_choose(n,m)*gsl_sf_fact(m);
      gsl_sf_result He;

      status = gsl_sf_hermite_prob_e(n - m, x, &He);
      if (status == GSL_SUCCESS)
        {
          result->val = He.val*f;
          result->err = He.err*f + GSL_DBL_EPSILON*fabs(result->val);
        }
      else
        {
          result->val = He.val;
          result->err = GSL_POSINF;
        }

      return status;
    }
}

double
gsl_sf_hermite_prob_deriv(const int m, const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_prob_deriv_e(m, n, x, &result));
}

/* evaluates the physicists' Hermite polynomial of order n at position x */
int
gsl_sf_hermite_e(const int n, const double x, gsl_sf_result * result)
{
  if (n < 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (n == 0)
    {
      result->val = 1.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
  else if (n == 1)
    {
      result->val = 2.0*x;
      result->err = 0.;
      return GSL_SUCCESS;
    }
  else if (x == 0.0)
    {
      if (GSL_IS_ODD(n))
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }
      else
        {
          /* for n even, H_n(0) = (-2)^{n/2} (n - 1)!! */

          int status = GSL_SUCCESS;
          int m = n >> 1;

          if (n - 1 > GSL_SF_DOUBLEFACT_NMAX)  /* test if (n - 1)!! will overflow */
            {
              status = GSL_EOVRFLW;
              result->val = GSL_IS_ODD(m) ? GSL_NEGINF : GSL_POSINF;
              result->err = GSL_POSINF;
            }
          else
            {
              double f = gsl_pow_int(2.0, m);

              gsl_sf_doublefact_e(n - 1, result);

              if (result->val > 0.9 * GSL_DBL_MAX / f) /* test if 2^{n/2} * (n-1)!! will overflow */
                {
                  status = GSL_EOVRFLW;
                  result->val = GSL_IS_ODD(m) ? GSL_NEGINF : GSL_POSINF;
                  result->err = GSL_POSINF;
                }
              else
                {
                  result->val *= f;
                  result->err *= f;
                  if (GSL_IS_ODD(m))
                    result->val = -result->val;
                }
            }

          return status;
        }
    }
  else
    {
      /* upward recurrence: H_{n+1} = 2x H_n - 2n H_{n-1} */

      int status = GSL_SUCCESS;
      const double two_x = 2.0 * x;
      const double abs_two_x = fabs(two_x);
      const double thresh1 = abs_two_x > 1.0 ? 0.9 * GSL_DBL_MAX / abs_two_x : GSL_DBL_MAX;
      const double thresh2 = 0.9 * GSL_DBL_MAX / 2.0;

      double p_n0 = 1.0;    /* H_0(x) */
      double p_n1 = two_x;  /* H_1(x) */
      double p_n = p_n1;

      double e_n0 = GSL_DBL_EPSILON;
      double e_n1 = 2.*fabs(x)*GSL_DBL_EPSILON;
      double e_n = e_n1;

      int j;

      for (j = 1; j <= n - 1; j++)
        {
          if (fabs(p_n1) > thresh1 ||    /* test if 2*x*p_n1 will overflow */
              fabs(p_n0) > thresh2 / j)  /* test if 2*j*p_n0 will overflow */
            {
              status = GSL_EOVRFLW;
              break;
            }

          p_n  = two_x*p_n1 - 2.0*j*p_n0;
          p_n0 = p_n1;
          p_n1 = p_n;

          e_n  = 2.*(fabs(x)*e_n1+j*e_n0);
          e_n0 = e_n1;
          e_n1 = e_n;
        }

      result->val = p_n;
      result->err = e_n + fabs(result->val)*GSL_DBL_EPSILON;

      return status;
    }
}

double
gsl_sf_hermite(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_e(n, x, &result));
}

/* Evaluates the m-th derivative of the physicists' Hermite polynomial of order n at position x.
 * The direct formula H^{(m)}_n = 2**m*n!/(n-m)!*H_{n-m}(x) (where H_j(x) is the j-th physicists' Hermite polynomial and H^{(m)}_j(x) its m-th derivative) is employed. */
int 
gsl_sf_hermite_deriv_e(const int m, const int n, const double x, gsl_sf_result * result)
{
  if (n < 0 || m < 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (n < m)
    {
      result->val = 0.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      double f = gsl_sf_choose(n,m)*gsl_sf_fact(m)*pow2(m);
      gsl_sf_result H;

      status = gsl_sf_hermite_e(n - m, x, &H);
      if (status == GSL_SUCCESS)
        {
          result->val = H.val*f;
          result->err = H.err*f + GSL_DBL_EPSILON*fabs(result->val);
        }
      else
        {
          result->val = H.val;
          result->err = GSL_POSINF;
        }

      return status;
    }
}

double
gsl_sf_hermite_deriv(const int m, const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_deriv_e(m, n, x, &result));
}

/* evaluates the Hermite function of order n at position x */
int
gsl_sf_hermite_func_e(const int n, const double x, gsl_sf_result * result)
{
  if (n < 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (x == 0.0)
    {
      if (GSL_IS_ODD(n))
        {
          result->val = 0.;
          result->err = 0.;
          return GSL_SUCCESS;
        }
      else
        {
          double f = (GSL_IS_ODD(n / 2) ? -1.0 : 1.0);
          int j;

          for(j = 1; j < n; j += 2)
            f *= sqrt(j / (j + 1.0));

          result->val = f / sqrt(M_SQRTPI);
          result->err = GSL_DBL_EPSILON * fabs(result->val);
          return GSL_SUCCESS;
        }
    }
  else if (n == 0)
    {
      result->val = exp(-0.5 * x * x) / sqrt(M_SQRTPI);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  else if (n == 1)
    {
      result->val = M_SQRT2 * x * exp(-0.5 * x * x) / sqrt(M_SQRTPI);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  else
    {
      /*
       * This algorithm is based on the modified recurrence algorithm
       * found in the appendix of:
       *
       * B. Bunck, BIT Numerical Mathematics, 49, 281 (2009)
       *
       * Numerical tests showed that this algorithm is more stable for
       * large x (x > 40) than the standard recurrence relation.
       * Accuracy is comparable to the recurrence relation method
       * for small x and all n.
       *
       * See:
       *
       * https://scicomp.stackexchange.com/questions/30896/generate-high-n-quantum-harmonic-oscillator-states-numerically
       *
       * for further discussion.
       */

      double hi2 = 1.0 / sqrt(M_SQRTPI); /* \hat{h}_0 */
      double hi1 = M_SQRT2 * x * hi2;    /* \hat{h}_1 */
      double hi = 0.0;
      double sum_log_scale = 0.0;
      double abshi;
      int i;

      for (i = 2; i <= n; ++i)
        {
          hi = sqrt(2.0 / i) * x * hi1 - sqrt((i - 1.0) / i) * hi2;
          hi2 = hi1;
          hi1 = hi;

          abshi = fabs(hi);
          if (abshi > 1.0)
            {
              double log_scale = RND(log(abshi));
              double scale = exp(-log_scale);

              hi *= scale;
              hi1 *= scale;
              hi2 *= scale;
              sum_log_scale += log_scale;
            }
        }

      result->val = hi * exp(-0.5 * x * x + sum_log_scale);
      result->err = n * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }
}

double
gsl_sf_hermite_func(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_func_e(n, x, &result));
}

/*
 * This algorithm is based on the contour integral algorithm of:
 *
 * B. Bunck, BIT Numerical Mathematics, 49, 281 (2009)
 *
 * It has O(sqrt(n)) complexity
 */

int
gsl_sf_hermite_func_fast_e(const int n, const double x, gsl_sf_result * result)
{
  if (n < 1000 || x == 0.0)
    {
      /* for small n, the recurrence method is faster and more accurate */
      return gsl_sf_hermite_func_e(n, x, result);
    }
  else
    {
      size_t j;
      const double k = sqrt(0.5*n);
      const size_t steps = (size_t) ceil(6.211 * sqrt(n));
      const double dt = M_PI/steps;
      const double invn2 = 1.0/(n*n);
      double ex, ex_e, cs, cs_e, sn, sn2, t;
      gsl_sf_result lngamma;
  
      if (n < 36)
        {
          gsl_sf_lnfact_e(n, &lngamma);
          lngamma.val *= 0.5;
          lngamma.err *= 0.5;
          t = 0.5*n*log(n) + 0.25*M_LNPI;
          cs = 0.5*n;
          lngamma.val += cs - t;
          lngamma.err += (cs + t)*GSL_DBL_EPSILON;
        }
      else
        {
          /* approximate ln(gamma_{n,k}) using Stirling's formula */
          lngamma.val = 0.25*log(2*n);
          lngamma.err = (lngamma.val  + ((((invn2/3360 + 1.0/2520)*invn2 + 1.0/720)*invn2) + 1.0/24)/n)*GSL_DBL_EPSILON;
          lngamma.val -= ((((invn2/3360 - 1.0/2520)*invn2 + 1.0/720)*invn2) - 1.0/24)/n;
        }
  
      ex = exp(lngamma.val - n - 0.5*x*x - 2*x*k);
      cs = (GSL_IS_ODD(n) ? -1 : 1);
      result->val = 0.5*ex*cs;
      result->err = 0.5*ex*(lngamma.err + (n + 0.5*x*x + fabs(2*x*k) + 1)*GSL_DBL_EPSILON);
      ex = exp(lngamma.val - n - 0.5*x*x + 2*x*k);
      result->val += 0.5*ex;
      result->err += 0.5*ex*(lngamma.err + (n + 0.5*x*x + fabs(2*x*k) + 1)*GSL_DBL_EPSILON);
  
      for (j = 1; j < steps; j++)
        {
          t = j*dt;
          cs = cos(t);
          ex = exp(lngamma.val - 0.5*x*x + (2*x*k - n*cs)*cs);
          ex_e = ex*(lngamma.err + GSL_DBL_EPSILON*(1 + 0.5*x*x + (fabs(2*x*k) + fabs(n*cs))*fabs(cs)));
          sn = sin(t);
          sn2 = sin(2*t);
          cs = cos(2*x*k*sn - 0.5*n*sn2 - n*t);
          cs_e = GSL_MIN(1.0+fabs(cs), GSL_DBL_EPSILON*(fabs(cs) + (fabs(2*x*k*sn) + fabs(0.5*n*sn2) + n*t)*fabs(sin(2*x*k*sn - 0.5*n*sn2 - n*t))));
          result->val += ex*cs;
          result->err += ex*cs_e + ex_e*fabs(cs) + GSL_DBL_EPSILON*fabs(ex*cs);
        }

      result->val *= M_1_PI*dt;
      result->err = M_1_PI*dt*result->err + GSL_DBL_EPSILON*fabs(result->val);

      return GSL_SUCCESS;
    }
}

double
gsl_sf_hermite_func_fast(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_func_fast_e(n, x, &result));
}

/* Evaluates all probabilists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_prob_array(const int nmax, const double x, double * result_array)
{
  if (nmax < 0)
    {
      GSL_ERROR ("domain error", GSL_EDOM);
    }
  else if (nmax == 0)
    {
      result_array[0] = 1.0;
      return GSL_SUCCESS;
    }
  else if (nmax == 1)
    {
      result_array[0] = 1.0;
      result_array[1] = x;
      return GSL_SUCCESS;
    }
  else
    {
      /* upward recurrence: He_{n+1} = x He_n - n He_{n-1} */

      int status = GSL_SUCCESS;
      const double abs_x = fabs(x);
      const double thresh1 = abs_x > 1.0 ? 0.9 * GSL_DBL_MAX / abs_x : GSL_DBL_MAX;
      const double thresh2 = 0.9 * GSL_DBL_MAX;

      double p_n0 = 1.0;    /* He_0(x) */
      double p_n1 = x;      /* He_1(x) */
      double p_n = p_n1;

      int j;

      result_array[0] = 1.0;
      result_array[1] = x;

      for (j = 1; j < nmax; j++)
        {
          if (fabs(p_n1) > thresh1 ||    /* test if x*p_n1 will overflow */
              fabs(p_n0) > thresh2 / j)  /* test if j*p_n0 will overflow */
            {
              status = GSL_EOVRFLW;
              break;
            }

          p_n  = x*p_n1 - j*p_n0;
          p_n0 = p_n1;
          p_n1 = p_n;

          result_array[j + 1] = p_n;
        }

      return status;
   }
}


/* Evaluates the m-th derivative of all probabilists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */

int
gsl_sf_hermite_prob_array_deriv(const int m, const int nmax, const double x, double * result_array)
{
  if (nmax < 0 || m < 0)
    {
      GSL_ERROR ("domain error", GSL_EDOM);
    }
  else if (m == 0)
    {
      gsl_sf_hermite_prob_array(nmax, x, result_array);
      return GSL_SUCCESS;
    }
  else if (nmax < m)
    {
      int j;

      for (j = 0; j <= nmax; j++)
        result_array[j] = 0.0;

      return GSL_SUCCESS;
    }
  else if (nmax == m)
    {
      int j;

      for (j = 0; j < m; j++)
        result_array[j] = 0.0;

      result_array[nmax] = gsl_sf_fact(m);
      return GSL_SUCCESS;
    }
  else if (nmax == m + 1)
    {
      int j;

      for (j = 0; j < m; j++)
        result_array[j] = 0.0;

      result_array[nmax-1] = gsl_sf_fact(m);
      result_array[nmax] = result_array[nmax-1]*(m+1)*x;
      return GSL_SUCCESS;
    }
  else
    {
      /* upward recurrence: He^{(m)}_{n+1} = (n+1)/(n-m+1)*(x He^{(m)}_n - n He^{(m)}_{n-1}) */

      double p_n0 = gsl_sf_fact(m);    /* He^{(m)}_{m}(x) */
      double p_n1 = p_n0*(m+1)*x;      /* He^{(m)}_{m+1}(x) */
      double p_n = p_n1;
      int j;

      for (j = 0; j < m; j++)
        result_array[j] = 0.0;

      result_array[m] = p_n0;
      result_array[m + 1] = p_n1;

      for (j = m + 1; j <= nmax - 1; j++)
        {
          p_n  = (x*p_n1 - j*p_n0) * (j + 1.0) / (j - m + 1.0);
          p_n0 = p_n1;
          p_n1 = p_n;
          result_array[j + 1] = p_n;
        }

      return GSL_SUCCESS;
    }
}

/* Evaluates all derivatives (starting from 0) up to the mmax-th derivative of the probabilists' Hermite polynomial of order n at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */

int
gsl_sf_hermite_prob_deriv_array(const int mmax, const int n, const double x, double * result_array)
{
  if (n < 0 || mmax < 0)
    {
      GSL_ERROR ("domain error", GSL_EDOM);
    }
  else if (n == 0)
    {
      int j;

      result_array[0] = 1.0;
      for (j = 1; j <= mmax; j++)
        result_array[j] = 0.0;

      return GSL_SUCCESS;
    }
  else if (n == 1 && mmax > 0)
    {
      int j;

      result_array[0] = x;
      result_array[1] = 1.0;

      for(j=2; j <= mmax; j++)
        result_array[j] = 0.0;

      return GSL_SUCCESS;
    }
  else if (mmax == 0)
    {
      result_array[0] = gsl_sf_hermite_prob(n,x);
      return GSL_SUCCESS;
    }
  else if (mmax == 1)
    {
      result_array[0] = gsl_sf_hermite_prob(n,x);
      result_array[1] = n*gsl_sf_hermite_prob(n-1,x);
      return GSL_SUCCESS;
    }
  else
    {
      /* upward recurrence */

      int k = GSL_MAX_INT(0, n - mmax);
      double p_n0 = gsl_sf_hermite_prob(k,x);        /* He_k(x) */
      double p_n1 = gsl_sf_hermite_prob(k+1,x);      /* He_{k+1}(x) */
      double p_n  = p_n1;
      int j;

      for(j=n+1; j <= mmax; j++)
        result_array[j] = 0.0;

      result_array[GSL_MIN_INT(n,mmax)] = p_n0;
      result_array[GSL_MIN_INT(n,mmax)-1] = p_n1;

      for (j = GSL_MIN_INT(mmax,n)-1; j > 0; j--)
        {
          k++;
          p_n  = x*p_n1-k*p_n0;
          p_n0 = p_n1;
          p_n1 = p_n;
          result_array[j - 1] = p_n;
        }

      p_n = 1.0;
      for (j = 1; j <= GSL_MIN_INT(n, mmax); j++)
        {
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
gsl_sf_hermite_array(const int nmax, const double x, double * result_array)
{
  if(nmax < 0)
    {
      GSL_ERROR ("domain error", GSL_EDOM);
    }
  else if (nmax == 0)
    {
      result_array[0] = 1.0;
      return GSL_SUCCESS;
    }
  else if (nmax == 1)
    {
      result_array[0] = 1.0;
      result_array[1] = 2.0*x;
      return GSL_SUCCESS;
    }
  else
    {
      /* upward recurrence: H_{n+1} = 2x H_n - 2n H_{n-1} */

      int status = GSL_SUCCESS;
      const double two_x = 2.0 * x;
      const double abs_two_x = fabs(two_x);
      const double thresh1 = abs_two_x > 1.0 ? 0.9 * GSL_DBL_MAX / abs_two_x : GSL_DBL_MAX;
      const double thresh2 = 0.9 * GSL_DBL_MAX / 2.0;

      double p_n0 = 1.0;    /* H_0(x) */
      double p_n1 = two_x;  /* H_1(x) */
      double p_n = p_n1;

      int j;

      result_array[0] = 1.0;
      result_array[1] = 2.0*x;

      for (j = 1; j < nmax; j++)
        {
          if (fabs(p_n1) > thresh1 ||    /* test if 2*x*p_n1 will overflow */
              fabs(p_n0) > thresh2 / j)  /* test if 2*j*p_n0 will overflow */
            {
              status = GSL_EOVRFLW;
            }

          p_n  = two_x*p_n1 - 2.0*j*p_n0;
          p_n0 = p_n1;
          p_n1 = p_n;

          result_array[j + 1] = p_n;
        }

      return status;
    }
}

/* Evaluates the m-th derivative of all physicists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_array_deriv(const int m, const int nmax, const double x, double * result_array)
{
  if (nmax < 0 || m < 0)
    {
      GSL_ERROR ("domain error", GSL_EDOM);
    }
  else if (m == 0)
    {
      gsl_sf_hermite_array(nmax, x, result_array);
      return GSL_SUCCESS;
    }
  else if (nmax < m)
    {
      int j;

      for(j = 0; j <= nmax; j++)
        result_array[j] = 0.0;

      return GSL_SUCCESS;
    }
  else if (nmax == m)
    {
      int j;

      for(j = 0; j < m; j++)
        result_array[j] = 0.0;

      result_array[nmax] = pow2(m)*gsl_sf_fact(m);

      return GSL_SUCCESS;
    }
  else if (nmax == m + 1)
    {
      int j;

      for(j = 0; j < m; j++)
        result_array[j] = 0.0;

      result_array[nmax-1] = pow2(m)*gsl_sf_fact(m);
      result_array[nmax] = result_array[nmax-1]*2*(m+1)*x;
      return GSL_SUCCESS;
    }
  else
    {
      /* upward recurrence: H^{(m)}_{n+1} = 2(n+1)/(n-m+1)*(x H^{(m)}_n - n H^{(m)}_{n-1}) */

      double p_n0 = pow2(m)*gsl_sf_fact(m);  /* H^{(m)}_{m}(x) */
      double p_n1 = p_n0*2*(m+1)*x;          /* H^{(m)}_{m+1}(x) */
      double p_n;
      int j;

      for(j = 0; j < m; j++)
        result_array[j] = 0.0;

      result_array[m] = p_n0;
      result_array[m+1] = p_n1;

      for (j = m + 1; j < nmax; ++j)
        {
          p_n  = (x*p_n1 - j*p_n0) * 2 * (j + 1.0) / (j - m + 1.0);
          p_n0 = p_n1;
          p_n1 = p_n;
          result_array[j + 1] = p_n;
        }

      return GSL_SUCCESS;
    }
}


/* Evaluates all derivatives (starting from 0) up to the mmax-th derivative of the physicists' Hermite polynomial of order n at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_deriv_array(const int mmax, const int n, const double x, double * result_array)
{
  if (n < 0 || mmax < 0)
    {
      GSL_ERROR ("domain error", GSL_EDOM);
    }
  else if (n == 0)
    {
      int j;

      result_array[0] = 1.0;
      for(j = 1; j <= mmax; j++)
        result_array[j] = 0.0;

      return GSL_SUCCESS;
    }
  else if (n == 1 && mmax > 0)
    {
      int j;

      result_array[0] = 2*x;
      result_array[1] = 1.0;
      for (j = 2; j <= mmax; j++)
        result_array[j] = 0.0;

      return GSL_SUCCESS;
    }
  else if (mmax == 0)
    {
      result_array[0] = gsl_sf_hermite(n,x);
      return GSL_SUCCESS;
    }
  else if (mmax == 1)
    {
      result_array[0] = gsl_sf_hermite(n,x);
      result_array[1] = 2*n*gsl_sf_hermite(n - 1, x);
      return GSL_SUCCESS;
    }
  else
    {
      /* upward recurrence */

      int k = GSL_MAX_INT(0, n - mmax);
      double p_n0 = gsl_sf_hermite(k, x);        /* H_k(x) */
      double p_n1 = gsl_sf_hermite(k + 1, x);    /* H_{k+1}(x) */
      double p_n  = p_n1;
      int j;

      for (j = n + 1; j <= mmax; j++)
        result_array[j] = 0.0;

      result_array[GSL_MIN_INT(n,mmax)] = p_n0;
      result_array[GSL_MIN_INT(n,mmax)-1] = p_n1;

      for (j = GSL_MIN_INT(mmax, n) - 1; j > 0; j--)
        {
          k++;
          p_n  = 2*x*p_n1 - 2*k*p_n0;
          p_n0 = p_n1;
          p_n1 = p_n;
          result_array[j - 1] = p_n;
        }

      p_n = 1.0;
      for (j = 1; j <= GSL_MIN_INT(n,mmax); j++)
        {
          p_n *= 2.0 * (n - j + 1.0);
          result_array[j] *= p_n;
        }

      return GSL_SUCCESS;
    }
}


/* Evaluates the series sum_{j=0}^n a_j*H_j(x) with H_j being the j-th physicists' Hermite polynomial.
 * For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to physicists' Hermite polynomials is used. */
int
gsl_sf_hermite_series_e(const int n, const double x, const double * a, gsl_sf_result * result)
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
gsl_sf_hermite_series(const int n, const double x, const double * a)
{
  EVAL_RESULT(gsl_sf_hermite_series_e(n, x, a, &result));
}


/* Evaluates all Hermite functions up to order nmax at position x. The results are stored in result_array.
 * Since all polynomial orders are needed, upward recurrence is employed. */
int
gsl_sf_hermite_func_array(const int nmax, const double x, double * result_array)
{
  if (nmax < 0)
    {
      GSL_ERROR ("domain error", GSL_EDOM);
    }
  else if (nmax == 0)
    {
      result_array[0] = exp(-0.5*x*x)/sqrt(M_SQRTPI);
      return GSL_SUCCESS;
    }
  else if (nmax == 1)
    {
      result_array[0] = exp(-0.5*x*x)/sqrt(M_SQRTPI);
      result_array[1] = result_array[0]*M_SQRT2*x;
      return GSL_SUCCESS;
    }
  else
    {
      /* upward recurrence: Psi_{n+1} = sqrt(2/(n+1))*x Psi_n - sqrt(n/(n+1)) Psi_{n-1} */

      const double arg = -0.5 * x * x;
      double hi2 = 1.0 / sqrt(M_SQRTPI);
      double hi1 = M_SQRT2 * x * hi2;
      double hi = 0.0;
      double sum_log_scale = 0.0;
      double abshi;
      int i;

      result_array[0] = exp(arg) * hi2;
      result_array[1] = result_array[0] * M_SQRT2 * x;

      for (i = 2; i <= nmax; ++i)
        {
          hi = sqrt(2.0 / i) * x * hi1 - sqrt((i - 1.0) / i) * hi2;
          hi2 = hi1;
          hi1 = hi;

          abshi = fabs(hi);
          if (abshi > 1.0)
            {
              double log_scale = RND(log(abshi));
              double scale = exp(-log_scale);

              hi *= scale;
              hi1 *= scale;
              hi2 *= scale;
              sum_log_scale += log_scale;
            }

          result_array[i] = hi * exp(arg + sum_log_scale);
        }

      return GSL_SUCCESS;
    }
}

/* Evaluates the series sum_{j=0}^n a_j*Psi_j(x) with Psi_j being the j-th Hermite function.
 * For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to Hermite functions is used. */

int
gsl_sf_hermite_func_series_e(const int n, const double x, const double * a, gsl_sf_result * result)
{
  if (n < 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (n == 0)
    {
      result->val = a[0]*exp(-0.5*x*x)/sqrt(M_SQRTPI);
      result->err = GSL_DBL_EPSILON*fabs(result->val);
      return GSL_SUCCESS;
    }
  else if (n == 1)
    {
      result->val = (a[0]+a[1]*M_SQRT2*x)*exp(-0.5*x*x)/sqrt(M_SQRTPI);
      result->err = 2.*GSL_DBL_EPSILON*(fabs(a[0])+fabs(a[1]*M_SQRT2*x))*exp(-0.5*x*x)/sqrt(M_SQRTPI);
      return GSL_SUCCESS;
    }
  else
    {
      /* downward recurrence: b_n = a_n + sqrt(2/(n+1))*x b_{n+1} - sqrt((n+1)/(n+2)) b_{n+2} */

      double b0 = 0.;
      double b1 = 0.;
      double btmp = 0.;

      double e0 = 0.;
      double e1 = 0.;
      double etmp = e1;

      int j;

      for (j = n; j >= 0; j--)
        {
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
  if(m < 0 || n < 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (m == 0)
    {
      return gsl_sf_hermite_func_e(n, x, result);
    }
  else if (m == 1)
    {
      double hi2 = 1.0 / sqrt(M_SQRTPI);
      double hi1 = M_SQRT2 * x * hi2;
      double hi = 0.0;
      double sum_log_scale = 0.0;
      double abshi;
      int i;

      for (i = 2; i <= n; ++i)
        {
          hi = sqrt(2.0 / i) * x * hi1 - sqrt((i - 1.0) / i) * hi2;
          hi2 = hi1;
          hi1 = hi;

          abshi = fabs(hi);
          if (abshi > 1.0)
            {
              double log_scale = RND(log(abshi));
              double scale = exp(-log_scale);

              hi *= scale;
              hi1 *= scale;
              hi2 *= scale;
              sum_log_scale += log_scale;
            }
        }

      /* psi'_n(x) = sqrt(2 n) psi_{n-1} - x psi_n */
      result->val = (sqrt(2.0*n) * hi2 - x * hi) * exp(-0.5 * x * x + sum_log_scale);
      result->err = n * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }
  else
    {
      int j;
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
        f *= sqrt(2.*j);

      if (m > n)
        {
          f = (GSL_IS_ODD(m-n)?-f:f);

          for (j=0;j<GSL_MIN_INT(n,m-n);j++)
            f *= (m-j)/(j+1.);
        }

      for (j=1;j<=m-n;j++)
        {
          b = x*h1-j*h0;
          h0 = h1;
          h1 = b;
	
          b  = (fabs(x)*eh1+j*eh0);
          eh0 = eh1;
          eh1 = b;
        }

      b = 0.;
      for (j=1;j<=n-m;j++)
        {
          b = (M_SQRT2*x*p1-sqrt(j)*p0)/sqrt(j+1.);
          p0 = p1;
          p1 = b;
	
          b  = (M_SQRT2*fabs(x)*ep1+sqrt(j)*ep0)/sqrt(j+1.);
          ep0 = ep1;
          ep1 = b;
        }

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
        }

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

/*
 * Computes the s-th zero the probabilists' Hermite polynomial of order n.
 * A Newton iteration using a continued fraction representation adapted from:
 *
 * [E.T. Whittaker (1914), On the continued fractions which represent the
 * functions of Hermite and other functions defined by differential equations,
 * Proceedings of the Edinburgh Mathematical Society, 32, 65-74]
 *
 * is performed with the initial approximation from
 *
 * [Arpad Elbert and Martin E. Muldoon, Approximations for zeros of Hermite
 * functions, pp. 117-126 in D. Dominici and R. S. Maier, eds, "Special Functions
 * and Orthogonal Polynomials", Contemporary Mathematics, vol 471 (2008)]
 *
 * refined via the bisection method.
 */

int
gsl_sf_hermite_prob_zero_e(const int n, const int s, gsl_sf_result * result)
{
  if (n <= 0 || s < 0 || s > n/2)
    {
      DOMAIN_ERROR(result);
    }
  else if (s == 0)
    {
      if (GSL_IS_ODD(n) == 1)
        {
          result->val = 0.;
          result->err = 0.;
          return GSL_SUCCESS;
        }
      else
        {
          DOMAIN_ERROR(result);
        }
    }
  else if (n == 2)
    {
      result->val = 1.;
      result->err = 0.;
      return GSL_SUCCESS;
    }
  else if (n < 21)
    {
      result->val = He_zero_tab[(GSL_IS_ODD(n)?n/2:0)+((n/2)*(n/2-1))+s-2];
      result->err = GSL_DBL_EPSILON*(result->val);
      return GSL_SUCCESS;
    }
  else
    {
      double d = 1., x = 1., x0 = 1.;
      int j;
      x = H_zero_init(n,s) * M_SQRT2;
      do
        {
          x0 = x;
          d = 0.;
          for (j=1; j<n; j++)
            d = j/(x-d);

          x -= (x-d)/n;

          /* gsl_fcmp can be used since the smallest zero approaches 1/sqrt(n) or 1/sqrt((n-1)/3.)
           * for large n and thus all zeros are non-zero (except for the trivial case handled above) */
        }
      while (gsl_fcmp(x, x0, 10*GSL_DBL_EPSILON) != 0);

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

/*
 * Computes the s-th zero the physicists' Hermite polynomial of order n, thus also
 * the s-th zero of the Hermite function of order n. A Newton iteration using a continued
 * fraction representation adapted from:
 *
 * [E.T. Whittaker (1914), On the continued fractions which represent the functions of Hermite
 * and other functions defined by differential equations, Proceedings of the Edinburgh Mathematical
 * Society, 32, 65-74]
 *
 * An initial approximation is used from:
 *
 * [Arpad Elbert and Martin E. Muldoon, Approximations for zeros of Hermite functions,
 * pp. 117-126 in D. Dominici and R. S. Maier, eds, "Special Functions and Orthogonal Polynomials",
 * Contemporary Mathematics, vol 471 (2008)]
 *
 * which is refined via the bisection method.
 */

int
gsl_sf_hermite_zero_e(const int n, const int s, gsl_sf_result * result)
{
  if (n <= 0 || s < 0 || s > n/2)
    {
      DOMAIN_ERROR(result);
    }
  else if (s == 0)
    {
      if (GSL_IS_ODD(n) == 1)
        {
          result->val = 0.;
          result->err = 0.;
          return GSL_SUCCESS;
        }
      else
        {
          DOMAIN_ERROR(result);
        }
    }
  else if (n == 2)
    {
      result->val = M_SQRT1_2;
      result->err = 0.;
      return GSL_SUCCESS;
    }
  else if (n < 21)
    {
      result->val = H_zero_tab[(GSL_IS_ODD(n)?n/2:0)+((n/2)*(n/2-1))+s-2];
      result->err = GSL_DBL_EPSILON*(result->val);
      return GSL_SUCCESS;
    }
  else
    {
      double d = 1., x = 1., x0 = 1.;
      int j;

      x = H_zero_init(n,s);
      do
        {
          x0 = x;
          d = 0.;

          for (j=1; j<n; j++)
            d = 2*j/(2.*x-d);

          x -= (2*x-d)*0.5/n;

        /* gsl_fcmp can be used since the smallest zero approaches 1/sqrt(n) or 1/sqrt((n-1)/3.)
         * for large n and thus all zeros are non-zero (except for the trivial case handled above) */
        }
      while (gsl_fcmp(x, x0, 10*GSL_DBL_EPSILON) != 0);

      result->val = x;
      result->err = 2*GSL_DBL_EPSILON*x + fabs(x-x0);

      return GSL_SUCCESS;
    }
}

double
gsl_sf_hermite_zero(const int n, const int s)
{
  EVAL_RESULT(gsl_sf_hermite_zero_e(n, s, &result));
}

int
gsl_sf_hermite_func_zero_e(const int n, const int s, gsl_sf_result * result)
{
  return gsl_sf_hermite_zero_e(n, s, result);
}

double
gsl_sf_hermite_func_zero(const int n, const int s)
{
  EVAL_RESULT(gsl_sf_hermite_func_zero_e(n, s, &result));
}

#ifndef GSL_DISABLE_DEPRECATED

int
gsl_sf_hermite_phys_e(const int n, const double x, gsl_sf_result * result)
{
  return gsl_sf_hermite_e(n, x, result);
}

double
gsl_sf_hermite_phys(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_phys_e(n, x, &result));
}

int 
gsl_sf_hermite_phys_der_e(const int m, const int n, const double x, gsl_sf_result * result)
{
  return gsl_sf_hermite_deriv_e(m, n, x, result);
}

double
gsl_sf_hermite_phys_der(const int m, const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_phys_der_e(m, n, x, &result));
}

int
gsl_sf_hermite_phys_array(const int nmax, const double x, double * result_array)
{
  return gsl_sf_hermite_array(nmax, x, result_array);
}

int
gsl_sf_hermite_phys_series_e(const int n, const double x, const double * a, gsl_sf_result * result)
{
  return gsl_sf_hermite_series_e(n, x, a, result);
}

double
gsl_sf_hermite_phys_series(const int n, const double x, const double * a)
{
  EVAL_RESULT(gsl_sf_hermite_phys_series_e(n, x, a, &result));
}

int
gsl_sf_hermite_phys_array_der(const int m, const int nmax, const double x, double * result_array)
{
  return gsl_sf_hermite_array_deriv(m, nmax, x, result_array);
}

int
gsl_sf_hermite_phys_der_array(const int mmax, const int n, const double x, double * result_array)
{
  return gsl_sf_hermite_deriv_array(mmax, n, x, result_array);
}

int
gsl_sf_hermite_phys_zero_e(const int n, const int s, gsl_sf_result * result)
{
  return gsl_sf_hermite_zero_e(n, s, result);
}

double
gsl_sf_hermite_phys_zero(const int n, const int s)
{
  EVAL_RESULT(gsl_sf_hermite_zero_e(n, s, &result));
}

int
gsl_sf_hermite_prob_array_der(const int m, const int nmax, const double x, double * result_array)
{
  return gsl_sf_hermite_prob_array_deriv(m, nmax, x, result_array);
}

int
gsl_sf_hermite_prob_der_array(const int mmax, const int n, const double x, double * result_array)
{
  return gsl_sf_hermite_prob_deriv_array(mmax, n, x, result_array);
}

int
gsl_sf_hermite_prob_der_e(const int m, const int n, const double x, gsl_sf_result * result)
{
  return gsl_sf_hermite_prob_deriv_e(m, n, x, result);
}

double
gsl_sf_hermite_prob_der(const int m, const int n, const double x)
{
  EVAL_RESULT(gsl_sf_hermite_prob_deriv_e(m, n, x, &result));
}

#endif /* !GSL_DISABLE_DEPRECATED */
