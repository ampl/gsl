/* specfunc/legendre_source.c
 * 
 * Copyright (C) 2009-2013 Patrick Alken
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

/* define various macros for functions below */

#define CONCAT2x(a,b)   a ## _ ## b 
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c

#if defined(LEGENDRE)
#define FUNCTION(dir,name) CONCAT2x(dir,name)
#define OUTPUT result_array
#define OUTPUT_ARG double result_array[]

#elif defined(LEGENDRE_DERIV)
#define FUNCTION(dir,name) CONCAT3x(dir,deriv,name)
#define OUTPUT result_array, result_deriv_array
#define OUTPUT_ARG double result_array[], double result_deriv_array[]

#elif defined(LEGENDRE_DERIV_ALT)
#define FUNCTION(dir,name) CONCAT3x(dir,deriv_alt,name)
#define OUTPUT result_array, result_deriv_array
#define OUTPUT_ARG double result_array[], double result_deriv_array[]
#define LEGENDRE_DERIV

#elif defined(LEGENDRE_DERIV2)
#define FUNCTION(dir,name) CONCAT3x(dir,deriv2,name)
#define OUTPUT result_array, result_deriv_array, result_deriv2_array
#define OUTPUT_ARG double result_array[], double result_deriv_array[], double result_deriv2_array[]
#define LEGENDRE_DERIV

#elif defined(LEGENDRE_DERIV2_ALT)
#define FUNCTION(dir,name) CONCAT3x(dir,deriv2_alt,name)
#define OUTPUT result_array, result_deriv_array, result_deriv2_array
#define OUTPUT_ARG double result_array[], double result_deriv_array[], double result_deriv2_array[]
#define LEGENDRE_DERIV
#define LEGENDRE_DERIV2
#define LEGENDRE_DERIV_ALT

#endif

static int FUNCTION (legendre, array_schmidt_e)
(const size_t lmax, const double x, const double csphase, OUTPUT_ARG);
static int FUNCTION(legendre, array_none_e)
(const size_t lmax, const double x, const double csphase, OUTPUT_ARG);

/*
gsl_sf_legendre_array()

Inputs: norm                - normlization type
        lmax                - maximum degree
        x                   - input argument
        result_array        - (output) normalized P_{lm}
        result_deriv_array  - (output) normalized P'_{lm}
        result_deriv2_array - (output) normalized P''_{lm}
*/

int
FUNCTION (gsl_sf_legendre, array)
(const gsl_sf_legendre_t norm, const size_t lmax, const double x,
 OUTPUT_ARG)
{
  int s = FUNCTION (gsl_sf_legendre, array_e)(norm, lmax, x, 1.0, OUTPUT);
  return s;
}

/*
gsl_sf_legendre_array_e()

Inputs: norm                - normlization type
        lmax                - maximum degree
        x                   - input argument
        csphase             - Condon-Shortley phase
        result_array        - (output) normalized P_{lm}
        result_deriv_array  - (output) normalized P'_{lm}
        result_deriv2_array - (output) normalized P''_{lm}
*/

int
FUNCTION (gsl_sf_legendre, array_e)
(const gsl_sf_legendre_t norm, const size_t lmax, const double x,
 const double csphase, OUTPUT_ARG)
{
  int s;
  const size_t nlm = gsl_sf_legendre_nlm(lmax);
#if !defined(LEGENDRE_DERIV_ALT)
  size_t i;
#if defined(LEGENDRE_DERIV) || defined(LEGENDRE_DERIV2)
  const double u = sqrt((1.0 - x) * (1.0 + x));
  const double uinv = 1.0 / u;
#endif
#if defined(LEGENDRE_DERIV2)
  const double uinv2 = uinv * uinv;
#endif
#endif
  double fac1 = 0.0, fac2 = 0.0; /* normalization factors */

  if (norm == GSL_SF_LEGENDRE_NONE)
    {
      /* compute P_{lm}(x) */
      s = FUNCTION(legendre,array_none_e)(lmax, x, csphase, OUTPUT);
    }
  else
    {
      /* compute S_{lm}(x) */
      s = FUNCTION(legendre,array_schmidt_e)(lmax, x, csphase, OUTPUT);
    }

#if !defined(LEGENDRE_DERIV_ALT)
  /* scale derivative arrays to recover P'(x) and P''(x) */
  for (i = 0; i < nlm; ++i)
    {
#if defined(LEGENDRE_DERIV2)
      double dp = result_deriv_array[i];
      double d2p = result_deriv2_array[i];

      result_deriv2_array[i] = (d2p - x * uinv * dp) * uinv2;
#endif
#if defined(LEGENDRE_DERIV)
      result_deriv_array[i] *= -uinv;
#endif
    }
#endif

  /* apply scaling for requested normalization */
  if (norm == GSL_SF_LEGENDRE_SCHMIDT || norm == GSL_SF_LEGENDRE_NONE)
    {
      return s;
    }
  else if (norm == GSL_SF_LEGENDRE_SPHARM)
    {
      fac1 = 1.0 / sqrt(4.0 * M_PI);
      fac2 = 1.0 / sqrt(8.0 * M_PI);
    }
  else if (norm == GSL_SF_LEGENDRE_FULL)
    {
      fac1 = 1.0 / sqrt(2.0);
      fac2 = 1.0 / sqrt(4.0);
    }

  /*
   * common code for different normalizations
   * P_{l0} = fac1 * sqrt(2l + 1) * S_{l0}
   * P_{lm} = fac2 * sqrt(2l + 1) * S_{lm}, m > 0
   */
  {
    size_t l, m;
    size_t twoellp1 = 1; /* 2l + 1 */
    double *sqrts = &(result_array[nlm]);

    for (l = 0; l <= lmax; ++l)
      {
        result_array[gsl_sf_legendre_array_index(l, 0)] *=
          sqrts[twoellp1] * fac1;
#if defined(LEGENDRE_DERIV)
        result_deriv_array[gsl_sf_legendre_array_index(l, 0)] *=
          sqrts[twoellp1] * fac1;
#endif
#if defined(LEGENDRE_DERIV2)
        result_deriv2_array[gsl_sf_legendre_array_index(l, 0)] *=
          sqrts[twoellp1] * fac1;
#endif

        for (m = 1; m <= l; ++m)
          {
            result_array[gsl_sf_legendre_array_index(l, m)] *=
              sqrts[twoellp1] * fac2;
#if defined(LEGENDRE_DERIV)
            result_deriv_array[gsl_sf_legendre_array_index(l, m)] *=
              sqrts[twoellp1] * fac2;
#endif
#if defined(LEGENDRE_DERIV2)
            result_deriv2_array[gsl_sf_legendre_array_index(l, m)] *=
              sqrts[twoellp1] * fac2;
#endif
          }

        twoellp1 += 2;
      }
  }

  return s;
}

/*
legendre,array_schmidt_e()
  This routine computes Schmidt semi-normalized associated
Legendre polynomials and their first and second derivatives.

Inputs: lmax                - maximum order
        x                   - legendre argument in [-1,1]
        csphase             - -1.0 to include CS phase (-1)^m,
                               1.0 to not include
        result_array        - (output) where to store P_{lm}(x) values
        result_deriv_array  - (output) where to store
                              d/dtheta P_{lm}(x) values
        result_deriv2_array - (output) where to store
                              d^2/dtheta^2 P_{lm}(x) values
*/

static int
FUNCTION(legendre, array_schmidt_e)
(const size_t lmax, const double x, const double csphase, OUTPUT_ARG)
{
  if (x > 1.0 || x < -1.0)
    {
      GSL_ERROR("x is outside [-1,1]", GSL_EDOM);
    }
#if defined(LEGENDRE_DERIV) || defined(LEGENDRE_DERIV2)
  else if (fabs(x) == 1.0)
    {
      GSL_ERROR("x cannot equal 1 or -1 for derivative computation", GSL_EDOM);
    }
#endif
  else if (csphase != 1.0 && csphase != -1.0)
    {
      GSL_ERROR("csphase has invalid value", GSL_EDOM);
    }
  else
    {
      const double eps = 1.0e-280;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
#if defined(LEGENDRE_DERIV)
      const double uinv = 1.0 / u;
#endif
#if defined(LEGENDRE_DERIV2)
      const double uinv2 = 1.0 / u / u;
#endif
#if defined(LEGENDRE_DERIV) || defined(LEGENDRE_DERIV2)
      const double xbyu = x * uinv; /* x / u */
#endif
      size_t l, m;
      size_t k, idxmm;
      double plm, /* eps * S(l,m) / u^m */
             pmm; /* eps * S(m,m) / u^m */
      double rescalem;
      double pm1, /* S(l-1,m) */
             pm2; /* S(l-2,m) */
      size_t nlm = gsl_sf_legendre_nlm(lmax);
      double *sqrts = &(result_array[nlm]);

      /* precompute square root factors for recurrence */
      legendre_sqrts(lmax, sqrts);

      /* initial values S(0,0) and S(1,0) */
      pm2 = 1.0; /* S(0,0) */
      pm1 = x;   /* S(1,0) */

      result_array[0] = pm2;
#if defined(LEGENDRE_DERIV)
      result_deriv_array[0] = 0.0;
#endif
#if defined(LEGENDRE_DERIV2)
      result_deriv2_array[0] = 0.0;
#endif

      if (lmax == 0)
        return GSL_SUCCESS;

      result_array[1] = pm1;
#if defined(LEGENDRE_DERIV)
      result_deriv_array[1] = -u;
#endif
#if defined(LEGENDRE_DERIV2)
      result_deriv2_array[1] = -x;
#endif

      /* Compute S(l,0) for l=2..lmax, no scaling required */

      k = 1; /* idx(1,0) */
      for (l = 2; l <= lmax; ++l)
        {
          double linv = 1.0 / (double)l;

          k += l;  /* idx(l,m) = idx(l-1,m) + l */

          plm = (2.0 - linv) * x * pm1 - (1.0 - linv) * pm2;
          result_array[k] = plm;
#if defined(LEGENDRE_DERIV)
          result_deriv_array[k] = uinv * l * (x * plm - pm1);
#endif
#if defined(LEGENDRE_DERIV2)
          result_deriv2_array[k] = -(double) l * (l + 1.0) * plm -
                                   xbyu * result_deriv_array[k];
#endif
          pm2 = pm1;
          pm1 = plm;
        }

      /* Compute S(m,m), S(m+1,m) and S(l,m) */

      /*
       * pi_m = Prod_{i=2}^m sqrt[ (2m - 1) / (2m) ]
       * but pi_1 = 1.0, so initialize to sqrt(2) so that
       * the first m = 1 iteration of the loop will reset it
       * to 1.0. Starting with m = 2 it will begin accumulating
       * the correct terms.
       *
       * pmm = S(m,m) * eps / u^m = pi_m
       */
      pmm = sqrt(2.0) * eps;

      rescalem = 1.0 / eps;
      idxmm = 0; /* tracks idx(m,m), initialize to idx(0,0) */

      for (m = 1; m < lmax; ++m)
        {
          /* rescalem = u^m / eps */
          rescalem *= u;

          /*
           * compute:
           * S(m,m) = u * sqrt((2m - 1) / (2m)) S(m-1,m-1) = u^m * pi_m
           * d_t S(m,m) = m * x / u * S(m,m)
           */

          idxmm += m + 1; /* idx(m,m) = idx(m-1,m-1) + m + 1 */
          pmm *= csphase * sqrts[2 * m - 1] / sqrts[2 * m]; /* S(m,m) * eps / u^m */
          result_array[idxmm] = pmm * rescalem;
#if defined(LEGENDRE_DERIV)
          result_deriv_array[idxmm] = m * xbyu * result_array[idxmm];
#endif
#if defined(LEGENDRE_DERIV2)
          result_deriv2_array[idxmm] =
            m * (uinv2 * m - (m + 1.0)) * result_array[idxmm] -
            xbyu * result_deriv_array[idxmm];
#endif
          pm2 = pmm;

          /*
           * compute:
           * S(m+1,m) = sqrt(2 * m + 1) * x * S(m,m)
           * d_t S(m+1,m) = 1/u * ((m+1)*x*S(m+1,m) - sqrt(2*m+1)*S(m,m))
           */

          k = idxmm + m + 1; /* idx(m+1,m) = idx(m,m) + m + 1 */
          pm1 = x * sqrts[2 * m + 1] * pm2;
          result_array[k] = pm1 * rescalem;
#if defined(LEGENDRE_DERIV)
          result_deriv_array[k] =
            uinv * ((m + 1.0) * x * result_array[k] -
                    sqrts[2 * m + 1] * result_array[idxmm]);
#endif
#if defined(LEGENDRE_DERIV2)
          result_deriv2_array[k] =
            (m * m * uinv2 - (m + 1.0) * (m + 2.0)) * result_array[k] -
            xbyu * result_deriv_array[k];
#endif

          /* compute S(l,m) for l=m+2...lmax */
          for (l = m + 2; l <= lmax; ++l)
            {
              k += l; /* idx(l,m) = idx(l-1,m) + l */
              plm =
                (2*l - 1) / sqrts[l + m] / sqrts[l - m] * x * pm1 -
                sqrts[l - m - 1] * sqrts[l + m - 1] /
                sqrts[l + m] / sqrts[l - m] * pm2;
              result_array[k] = plm * rescalem;
#if defined(LEGENDRE_DERIV)
              result_deriv_array[k] =
                uinv * (l * x * result_array[k] -
                        sqrts[l + m] * sqrts[l - m] * result_array[k - l]);
#endif
#if defined(LEGENDRE_DERIV2)
              result_deriv2_array[k] =
                (m * m * uinv2 - l * (l + 1.0)) * result_array[k] -
                xbyu * result_deriv_array[k];
#endif
              pm2 = pm1;
              pm1 = plm;
            }
        } /* for (m = 1; m < lmax; ++m) */

      /* compute S(lmax,lmax) */

      rescalem *= u;
      idxmm += m + 1; /* idx(lmax,lmax) */
      pmm *= csphase * sqrts[2 * lmax - 1] / sqrts[2 * lmax];
      result_array[idxmm] = pmm * rescalem;
#if defined(LEGENDRE_DERIV)
      result_deriv_array[idxmm] = lmax * xbyu * result_array[idxmm];
#endif
#if defined(LEGENDRE_DERIV2)
      result_deriv2_array[idxmm] =
        lmax * (uinv2 * lmax - (lmax + 1.0)) * result_array[idxmm] -
        xbyu * result_deriv_array[idxmm];
#endif

      return GSL_SUCCESS;
    }
}

/*
legendre_array_none_e()
  This routine computes unnormalized associated Legendre polynomials
and their derivatives.

Inputs: lmax                - maximum order
        x                   - legendre argument in [-1,1]
        csphase             - -1.0 to include CS phase (-1)^m,
                               1.0 to not include
        result_array        - (output) where to store P_{lm}(x) values
        result_deriv_array  - (output) where to store
                              d/dtheta P_{lm}(x) values
        result_deriv2_array - (output) where to store
                              d^2/dtheta^2 P_{lm}(x) values
*/

static int
FUNCTION(legendre, array_none_e)
(const size_t lmax, const double x, const double csphase, OUTPUT_ARG)
{
  if (x > 1.0 || x < -1.0)
    {
      GSL_ERROR("x is outside [-1,1]", GSL_EDOM);
    }
#if defined(LEGENDRE_DERIV) || defined(LEGENDRE_DERIV2)
  else if (fabs(x) == 1.0)
    {
      GSL_ERROR("x cannot equal 1 or -1 for derivative computation", GSL_EDOM);
    }
#endif
  else if (csphase != 1.0 && csphase != -1.0)
    {
      GSL_ERROR("csphase has invalid value", GSL_EDOM);
    }
  else
    {
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
#if defined(LEGENDRE_DERIV)
      const double uinv = 1.0 / u;
#endif
#if defined(LEGENDRE_DERIV2)
      const double uinv2 = 1.0 / u / u;
#endif
#if defined(LEGENDRE_DERIV) || defined(LEGENDRE_DERIV2)
      const double xbyu = x * uinv; /* x / u */
#endif
      size_t l, m;
      size_t k, idxmm;
      double plm, pmm;
      double pm1,    /* P(l-1,m) */
             pm2;    /* P(l-2,m) */
      double twomm1; /* 2*m - 1 */

      /* initial values P(0,0) and P(1,0) */

      pm2 = 1.0; /* P(0,0) */
      pm1 = x;   /* P(1,0) */

      result_array[0] = pm2;
#if defined(LEGENDRE_DERIV)
      result_deriv_array[0] = 0.0;
#endif
#if defined(LEGENDRE_DERIV2)
      result_deriv2_array[0] = 0.0;
#endif

      if (lmax == 0)
        return 0;

      result_array[1] = pm1;
#if defined(LEGENDRE_DERIV)
      result_deriv_array[1] = -u;
#endif
#if defined(LEGENDRE_DERIV2)
      result_deriv2_array[1] = -x;
#endif

      /* Compute P(l,0) */

      k = 1;
      for (l = 2; l <= lmax; ++l)
        {
          k += l;
          plm = ((2*l - 1) * x * pm1 - (l - 1) * pm2) / (double) l;
          result_array[k] = plm;
#if defined(LEGENDRE_DERIV)
          result_deriv_array[k] = -(double)l * (pm1 - x * plm) * uinv;
#endif
#if defined(LEGENDRE_DERIV2)
          result_deriv2_array[k] = -(double) l * (l + 1.0) * plm -
                                   xbyu * result_deriv_array[k];
#endif
          pm2 = pm1;
          pm1 = plm;
        }

      /* Compute P(m,m), P(m+1,m) and P(l,m) */

      pmm = 1.0;
      twomm1 = -1.0; /* 2 * m - 1 */

      idxmm = 0; /* tracks idx(m,m), initialize to idx(0,0) */

      for (m = 1; m <= lmax - 1; ++m)
        {
          /*
           * compute
           *
           * P(m,m) = u * (2m - 1) P(m-1,m-1)
           * and
           * dP(m,m)/dtheta = m * x * P(m,m) / u
           */
          idxmm += m + 1;
          twomm1 += 2.0;
          pmm *= csphase * u * twomm1;
          result_array[idxmm] = pmm;
#if defined(LEGENDRE_DERIV)
          result_deriv_array[idxmm] = m * xbyu * pmm;
#endif
#if defined(LEGENDRE_DERIV2)
          result_deriv2_array[idxmm] =
            m * (uinv2 * m - (m + 1.0)) * result_array[idxmm] -
            xbyu * result_deriv_array[idxmm];
#endif
          pm2 = pmm;

          /*
           * compute
           *
           * P(m+1,m) = (2 * m + 1) * x * P(m,m)
           * and
           * dP(m+1,m)/dt = -[(2*m + 1) * P(m,m) - (m+1) * x * P(m+1,m)]/u
           */
          k = idxmm + m + 1;
          pm1 = x * pmm * (2*m + 1);
          result_array[k] = pm1;
#if defined(LEGENDRE_DERIV)
          result_deriv_array[k] = -uinv * ((2*m + 1) * pmm - (m + 1) * x * pm1);
#endif
#if defined(LEGENDRE_DERIV2)
          result_deriv2_array[k] =
            (m * m * uinv2 - (m + 1.0) * (m + 2.0)) * result_array[k] -
            xbyu * result_deriv_array[k];
#endif

          /* compute P(l,m) */
          for (l = m + 2; l <= lmax; ++l)
            {
              k += l;
              plm = ((2*l - 1) * x * pm1 - (l + m - 1) * pm2) /
                    (double) (l - m);
              result_array[k] = plm;
#if defined(LEGENDRE_DERIV)
              result_deriv_array[k] = -uinv * ((l + m) * pm1 - l * x * plm);
#endif
#if defined(LEGENDRE_DERIV2)
              result_deriv2_array[k] =
                (m * m * uinv2 - l * (l + 1.0)) * result_array[k] -
                xbyu * result_deriv_array[k];
#endif
              pm2 = pm1;
              pm1 = plm;
            }
        } /* for (m = 1; m <= lmax - 1; ++m) */

      /* compute P(lmax,lmax) */

      idxmm += m + 1;
      twomm1 += 2.0;
      pmm *= csphase * u * twomm1;
      result_array[idxmm] = pmm;
#if defined(LEGENDRE_DERIV)
      result_deriv_array[idxmm] = lmax * x * pmm * uinv;
#endif
#if defined(LEGENDRE_DERIV2)
      result_deriv2_array[idxmm] =
        lmax * (uinv2 * lmax - (lmax + 1.0)) * result_array[idxmm] -
        xbyu * result_deriv_array[idxmm];
#endif

      return GSL_SUCCESS;
    }
} /* legendre_array_none_e() */

#undef FUNCTION
#undef CONCAT2x
#undef CONCAT3x
#undef OUTPUT
#undef OUTPUT_ARG
#undef LEGENDRE_DERIV
#undef LEGENDRE_DERIV2
#undef LEGENDRE_DERIV_ALT
