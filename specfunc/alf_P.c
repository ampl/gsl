/* specfunc/alf_P.c
 * 
 * Copyright (C) 2009-2023 Patrick Alken
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
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_alf.h>

/*
 * The routines in this module compute associated Legendre functions
 * (ALFs) up to order and degree 2700, using the method described
 * in
 *
 * [1] S. A. Holmes and W. E. Featherstone, A unified approach
 *     to the Clenshaw summation and the recursive computation of very
 *     high degree and order normalised associated Legendre functions,
 *     Journal of Geodesy, 76, pg. 279-299, 2002.
 *
 * Computations of the alternative derivatives, d^k/dtheta^k Plm
 * use the methods described in,
 *
 * [2] W. Bosch, On the computation of derivatives of Legendre functions,
 *     Phys. Chem. Earth (A), 25, 9-11, pg. 655-659, 2000.
 *
 * Further information on ALFs can be found in
 *
 * [3] Abramowitz and Stegun, Handbook of Mathematical Functions,
 *     Chapter 8, 1972.
 */

/*
gsl_sf_alf_precompute()
  Precompute recurrence factors for ALFs

Inputs: norm         - ALF normalization
        lmax         - maximum degree
        mmax         - maximum order, must be <= lmax
        flags        - GSL_SF_ALF_FLG_xxx
        output_array - (output) output array, length determined by
                       gsl_sf_alf_array_size()

Notes:
1) Precomputed factors are stored at the back of output_array, skipping
nlm entries which will be used for the ALF values

2) The precomputed factors are:
alm:   length nlm
blm:   length nlm
cl:    length lmax+1
dl:    length lmax+1
sqrts: length 2*lmax+2

Special storage:

a(0,0) = alm[2*idx(0,0) = 0]           = P(0,0)
a(1,0) = alm[2*idx(1,0) = 2]           = unused
a(1,1) = alm[2*idx(1,1) = 2*(L+1)]     = unused
b(0,0) = alm[2*idx(0,0)+1 = 1]         = unused
b(1,0) = alm[2*idx(1,0)+1 = 3]         = unused
b(1,1) = alm[2*idx(1,1)+1 = 2*(L+1)+1] = unused
dl[0]                                  = csfac
*/

int
gsl_sf_alf_precompute(const gsl_sf_alf_t norm, const size_t lmax,
                      const size_t mmax, const size_t flags, double output_array[])
{
  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else
    {
      const double csfac = (flags & GSL_SF_ALF_FLG_CSPHASE) ? -1.0 : 1.0;
      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      double * alm = &output_array[nlm]; /* save nlm entries for the ALFs */
      double * cl = alm + 2 * nlm;
      double * dl = cl + lmax + 1;
      double * sqrts = dl + lmax + 1;
      size_t l, m, k;

      /* compute square root factors */
      for (l = 0; l <= 2 * lmax + 1; ++l)
        sqrts[l] = sqrt((double) l);

      /* store csfac in dl[0] which is not needed in the recurrence relations */
      dl[0] = csfac;

      if (norm == GSL_SF_ALF_SCHMIDT)
        {
          alm[0] = 1.0; /* S(0,0) */
          cl[0] = 1.0;

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = M_SQRT3;
          dl[1] = csfac;

          /* m = 0 terms */
          k = 2; /* idx(2,0) */
          for (l = 2; l <= lmax; ++l)
            {
              /* a_l^0 */
              alm[2*k] = 2.0 - 1.0 / (double) l;

              /* b_l^0 */
              alm[2*k + 1] = -(1.0 - 1.0 / (double) l);

              ++k;
            }

          /* m > 0 terms */
          k = lmax + 1; /* idx(1,1) */
          for (m = 1; m <= mmax; ++m)
            {
              /* alm and blm are unused for l=m and l=m+1 */
              k += 2;

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = ((2.0*l - 1.0) / sqrts[l + m]) / sqrts[l - m];

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrts[l + m - 1] / sqrts[l + m]) *
                                  (sqrts[l - m - 1] / sqrts[l - m]);

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrts[2 * l + 1];
              dl[l] = csfac * sqrt(1.0 - 0.5 / l);
            }
        }
      else if (norm == GSL_SF_ALF_FOURPI)
        {
          alm[0] = 1.0;    /* R(0,0) */
          cl[0] = M_SQRT3;

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = 2.236067977499789696409174; /* sqrt(5) */
          dl[1] = csfac * M_SQRT3;

          k = 2; /* idx(2,0) */
          for (m = 0; m <= mmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (sqrts[2 * l + 1] / sqrts[l + m]) *
                             (sqrts[2 * l - 1] / sqrts[l - m]);

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrts[l + m - 1] / sqrts[l + m]) *
                                  (sqrts[l - m - 1] / sqrts[l - m]) *
                                  (sqrts[2*l + 1] / sqrts[2*l - 3]);

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 3.0);
              dl[l] = csfac * sqrt(1.0 + 0.5 / l);
            }
        }
      else if (norm == GSL_SF_ALF_SPHARM)
        {
          alm[0] = 0.5 / M_SQRTPI;    /* Y(0,0) */
          cl[0] = M_SQRT3;

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = sqrt(5.0);
          dl[1] = csfac * (M_SQRT3 / M_SQRT2);

          k = 2; /* idx(2,0) */
          for (m = 0; m <= mmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (sqrts[2 * l + 1] / sqrts[l + m]) *
                             (sqrts[2 * l - 1] / sqrts[l - m]);

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrts[l + m - 1] / sqrts[l + m]) *
                                  (sqrts[l - m - 1] / sqrts[l - m]) *
                                  (sqrts[2*l + 1] / sqrts[2*l - 3]);

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 3.0);
              dl[l] = csfac * sqrt(1.0 + 0.5 / l);
            }
        }
      else if (norm == GSL_SF_ALF_FULL)
        {
          cl[0] = M_SQRT3;
          alm[0] = M_SQRT1_2;         /* N(0,0) */

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = sqrt(5.0);
          dl[1] = csfac * (M_SQRT3 / M_SQRT2);

          k = 2; /* idx(2,0) */
          for (m = 0; m <= mmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (sqrts[2 * l + 1] / sqrts[l + m]) *
                             (sqrts[2 * l - 1] / sqrts[l - m]);

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrts[l + m - 1] / sqrts[l + m]) *
                                  (sqrts[l - m - 1] / sqrts[l - m]) *
                                  (sqrts[2*l + 1] / sqrts[2*l - 3]);

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 3.0);
              dl[l] = csfac * sqrt(1.0 + 0.5 / l);
            }
        }
      else if (norm == GSL_SF_ALF_NONE)
        {
          cl[0] = 1.0;
          alm[0] = 1.0;         /* P(0,0) */

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = 3.0;
          dl[1] = csfac;

          k = 2; /* idx(2,0) */
          for (m = 0; m <= mmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (2.0 * l - 1.0) / ((double) (l - m));

                  /* b_l^m */
                  alm[2*k + 1] = -(l + m - 1.0) / ((double) (l - m));

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = 2.0 * l + 1.0;
              dl[l] = csfac * (2.0 * l - 1.0);
            }
        }
      else
        {
          GSL_ERROR ("unknown normalization", GSL_EDOM);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_array_size()
  Compute total size of array needed for array functions

size = nlm +     // for ALFs
       nlm +     // for alm factors
       nlm +     // for blm factors
       lmax+1    // for cl factors
       lmax+1    // for dl factors
       2*lmax+2  // for sqrt factors

Inputs: lmax - maximum degree
*/

size_t
gsl_sf_alf_array_size(const size_t lmax, const size_t mmax)
{
  const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
  const size_t size = 3 * nlm + 4 * (lmax + 1);
  return size;
}

/*
gsl_sf_alf_array()
  Compute array of associated Legendre functions at a
given point x

Inputs: lmax         - maximum degree
        mmax         - maximum order (<= lmax)
        x            - input point
        result_array - (input/output) output array of ALFs

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) result_array is indexed by gsl_sf_alf_array_index(l,m,lmax),
result_array[index(l,m,lmax)] = Plm(x)
*/

int
gsl_sf_alf_array(const size_t lmax, const size_t mmax, const double x,
                 double result_array[])
{
  const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
  const double * alm = &result_array[nlm];
  const double * cl = alm + 2 * nlm;
  const double * dl = cl + lmax + 1;

  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else
    {
      /* interior point */

      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const size_t Lp2 = lmax + 2;
      size_t l, m;
      double plm,      /* eps * P(l,m) */
             pmm;      /* eps * P(m,m) */
      double plmp1;    /* P(l+1,m) */
      size_t idxmm;    /* idx(m,m) */
      const double *al;

      /* initial values P(0,0) and P(1,0) */

      plm = alm[0];            /* P(0,0) */
      result_array[0] = plm;

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      plmp1 = cl[0] * x * plm; /* P(1,0) */
      result_array[1] = plmp1;

      /* Compute P(l,0), l=2:lmax */

      l = 2;        /* idx(2,0) */
      al = &alm[4]; /* a(2,0) */
      while (l < lmax)
        {
          plm   = (al[0] * x) * plmp1 + al[1] * plm;
          plmp1 = (al[2] * x) * plm   + al[3] * plmp1;

          result_array[l] = plm; result_array[l+1] = plmp1;
          l += 2; al += 4;
        }

      if (l == lmax)
        result_array[l] = (al[0] * x) * plmp1 + al[1] * plm;

      if (mmax == 0)
        return GSL_SUCCESS;

      /* compute P(m,m), P(m+1,m) and P(l,m) */

      pmm = result_array[0]; /* P(0,0) */
      idxmm = 0;             /* idx(0,0) */

      for (m = 1; m <= mmax; ++m)
        {
          double *ptr;

          /* compute P(m,m) = d_m * u * P(m-1,m-1) */
          idxmm += Lp2 - m; /* idx(m,m) = idx(m-1,m-1) + L + 2 - m */
          pmm *= dl[m] * u;
          result_array[idxmm] = pmm;
          plm = pmm;

          if (m + 1 <= lmax)
            {
              /* compute P(m+1,m) = c_m * x * P(m,m) */
              plmp1 = (cl[m] * x) * plm;
              result_array[idxmm + 1] = plmp1;

              /* compute P(l,m) for l=m+2:lmax */
              l = m + 2;
              ptr = result_array + idxmm - m; /* shift array pointer */
              al = &alm[2 * (idxmm + 2)];     /* a(m+2,m) */
              while (l < lmax)
                {
                  plm   = (al[0] * x) * plmp1 + al[1] * plm;
                  plmp1 = (al[2] * x) * plm   + al[3] * plmp1;

                  ptr[l]     = plm;
                  ptr[l + 1] = plmp1;

                  l += 2; al += 4;
                }

              if (l == lmax)
                {
                  plm = (al[0] * x) * plmp1 + al[1] * plm;
                  ptr[l] = plm;
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_deriv_array()
  Compute array of associated Legendre functions at a
given point x

Inputs: lmax               - maximum degree
        mmax               - maximum order (<= lmax)
        x                  - input point in (-1,1)
        result_array       - (input/output) output array of ALFs
        result_deriv_array - (input/output) output array of dPlm(x)/dx

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) output arrays are indexed by gsl_sf_alf_array_index(l,m,lmax),
result_array[index(l,m,lmax)] = Plm(x)
result_deriv_array[index(l,m,lmax)] = dPlm(x)/dx
*/

int
gsl_sf_alf_deriv_array(const size_t lmax, const size_t mmax, const double x,
                       double result_array[], double result_deriv_array[])
{
  const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
  const double * alm = &result_array[nlm];
  const double * cl = alm + 2 * nlm;
  const double * dl = cl + lmax + 1;

  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else if (x == -1.0 || x == 1.0) /* endpoints */
    {
      GSL_ERROR ("x cannot equal endpoints", GSL_EDOM);
    }
  else
    {
      /* interior point */

      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double ratio = x / u;
      const size_t Lp2 = lmax + 2;
      size_t l, m, k;
      double plm,      /* P(l,m) */
             pmm;      /* P(m,m) */
      double plmp1;    /* P(l+1,m) */
      double dplm,
             dplmp1,
             dpmm;
      size_t idxmm;    /* idx(m,m) */
      const double *al;

      /* initial values P(0,0) and dP(0,0)/dx */
      plm = alm[0];                /* P(0,0) */
      dplm = 0.0;                  /* dP(0,0)/dx */
      result_array[0] = plm;
      result_deriv_array[0] = dplm;

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      /* compute P(1,0) and dP(1,0)/dx */
      plmp1 = cl[0] * x * plm; /* P(1,0) */
      dplmp1 = cl[0] * plm;    /* dP(1,0)/dx */
      result_array[1] = plmp1;
      result_deriv_array[1] = dplmp1;

      /* compute P(l,0) and dP(l,0)/dx, l=2:lmax */

      k = 2; /* idx(2,0) */
      al = &alm[k << 1];
      for (l = 2; l < lmax; l += 2)
        {
          plm    = (al[0] * x) * plmp1 + al[1] * plm;
          dplm   = al[0] * (x * dplmp1 + plmp1) + al[1] * dplm;

          plmp1  = (al[2] * x) * plm + al[3] * plmp1;
          dplmp1 = al[2] * (x * dplm + plm) + al[3] * dplmp1;

          result_array[k] = plm;        result_array[k+1] = plmp1;
          result_deriv_array[k] = dplm; result_deriv_array[k+1] = dplmp1;
          k += 2; al += 4;
        }

      if (l == lmax)
        {
          result_array[k] = (al[0] * x) * plmp1 + al[1] * plm;
          result_deriv_array[k] = al[0] * (x * dplmp1 + plmp1) + al[1] * dplm;
        }

      if (mmax == 0)
        return GSL_SUCCESS;

      /* compute P(m,m), P(m+1,m) and P(l,m) */

      pmm = result_array[0]; /* P(0,0) */
      dpmm = 0.0;            /* dP(0,0)/dx */
      idxmm = 0;             /* idx(0,0) */

      for (m = 1; m <= mmax; ++m)
        {
          /* compute P(m,m) = d_m * u * P(m-1,m-1) */
          idxmm += Lp2 - m; /* idx(m,m) = idx(m-1,m-1) + L + 2 - m */
          dpmm = dl[m] * (-ratio * pmm + u * dpmm);
          pmm *= dl[m] * u;
          result_array[idxmm] = pmm;
          result_deriv_array[idxmm] = dpmm;
          plm = pmm;
          dplm = dpmm;

          if (m + 1 <= lmax)
            {
              /* compute P(m+1,m) = c_m * x * P(m,m) */
              k = idxmm + 1; /* idx(m+1,m) = idx(m,m) + 1 */
              plmp1 = (cl[m] * x) * plm;
              dplmp1 = cl[m] * (pmm + x * dpmm);
              result_array[k] = plmp1;
              result_deriv_array[k] = dplmp1;

              /* compute P(l,m) for l=m+2:lmax */
              ++k;
              al = &alm[k << 1];
              for (l = m + 2; l < lmax; l += 2)
                {
                  plm    = (al[0] * x) * plmp1 + al[1] * plm;
                  dplm   = al[0] * (x * dplmp1 + plmp1) + al[1] * dplm;

                  plmp1  = (al[2] * x) * plm + al[3] * plmp1;
                  dplmp1 = al[2] * (x * dplm + plm) + al[3] * dplmp1;

                  result_array[k] = plm; result_array[k + 1] = plmp1;
                  result_deriv_array[k] = dplm; result_deriv_array[k + 1] = dplmp1;

                  k += 2; al += 4;
                }

              if (l == lmax)
                {
                  result_array[k] = (al[0] * x) * plmp1 + al[1] * plm;
                  result_deriv_array[k] = al[0] * (x * dplmp1 + plmp1) + al[1] * dplm;
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_vsh_array()
  Compute array of associated Legendre functions and their
theta derivatives at a given point x suitable for vector spherical
harmonic expansions

Inputs: lmax               - maximum degree
        mmax               - maximum order
        x                  - input point
        result_array       - (input/output) output array of ALFs
                             Plm(x)            for m = 0
                             Plm(x)/sin(theta) for m > 0
        result_deriv_array - (output) output array of ALF derivatives
                             d/dtheta Plm(x)

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) result_array and result_deriv_array are indexed by gsl_sf_alf_array_index(l,m,lmax),
result_array[index(l,0,lmax)]       = Pl0(x)
result_array[index(l,m,lmax)]       = Plm(x)/sin(theta), m > 0
result_deriv_array[index(l,m,lmax)] = dPlm(x)/dtheta
*/

int
gsl_sf_alf_vsh_array(const size_t lmax, const size_t mmax, const double x,
                     double result_array[], double result_deriv_array[])
{
  const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
  const double * alm = &result_array[nlm];
  const double * cl = alm + 2 * nlm;
  const double * dl = cl + lmax + 1;

  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else
    {
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double usq = u * u;
      const size_t Lp1 = lmax + 1;
      size_t l, m;
      double plm,      /* P(l,m) */
             pmm;      /* P(m,m) */
      double plmp1;    /* P(l+1,m) */
      double dplm,
             dplmp1,
             dpmm;
      size_t idxmm;    /* idx(m,m) */
      const double *al;

      /* initial values P(0,0) and dP(0,0)/dtheta */
      plm = alm[0];                /* P(0,0) */
      dplm = 0.0;                  /* dP(0,0)/dtheta */
      result_array[0] = plm;
      result_deriv_array[0] = dplm;

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      /* compute P(1,0) and dP(1,0)/dtheta */
      plmp1 = cl[0] * x * plm;    /* P(1,0) */
      dplmp1 = -u * cl[0] * plm;  /* dP(1,0)/dtheta */
      result_array[1] = plmp1;
      result_deriv_array[1] = dplmp1;

      /* Compute P(l,0), l=2:lmax */

      l = 2;        /* idx(2,0) */
      al = &alm[4]; /* a(2,0) */
      while (l < lmax)
        {
          plm    = (al[0] * x) * plmp1 + al[1] * plm;
          dplm   = al[0] * (x * dplmp1 - u * plmp1) + al[1] * dplm;

          plmp1  = (al[2] * x) * plm + al[3] * plmp1;
          dplmp1 = al[2] * (x * dplm - u * plm) + al[3] * dplmp1;

          result_array[l] = plm; result_array[l+1] = plmp1;
          result_deriv_array[l] = dplm; result_deriv_array[l+1] = dplmp1;
          l += 2; al += 4;
        }

      if (l == lmax)
        {
          result_array[l] = (al[0] * x) * plmp1 + al[1] * plm;
          result_deriv_array[l] = al[0] * (x * dplmp1 - u * plmp1) + al[1] * dplm;
        }

      if (mmax == 0)
        return GSL_SUCCESS;

      /* compute P(1,1)/u and dP(1,1)/dtheta */

      pmm = dl[1] * result_array[0];       /* P(1,1) / sin(theta) */
      dpmm = dl[1] * x * result_array[0];  /* dP(1,1)/dtheta */
      idxmm = Lp1;                         /* idx(1,1) */
      result_array[idxmm] = pmm;
      result_deriv_array[idxmm] = dpmm;

      if (lmax == 1)
        return GSL_SUCCESS;

      /* compute P(m,m), P(m+1,m) and P(l,m) */

      for (m = 1; m <= mmax; ++m)
        {
          double *ptr, *dptr;

          plm = pmm;   /* P(m,m)/sin(theta) */
          dplm = dpmm; /* dP(m,m)/dtheta */

          if (m + 1 <= lmax)
            {
              /* compute P(m+1,m)/sint = c_m * x * P(m,m)/sint */
              plmp1 = (cl[m] * x) * plm;

              /* compute dP(m+1,m)/dtheta = c_m * (-u^2 P(m,m)/u + x dP(m,m)/dtheta) */
              dplmp1 = cl[m] * (-usq * plm + x * dplm);
              result_array[idxmm + 1] = plmp1;
              result_deriv_array[idxmm + 1] = dplmp1;
            }

          /* compute P(l,m) for l=m+2:lmax */
          l = m + 2;
          ptr = result_array + idxmm - m;        /* shift array pointers */
          dptr = result_deriv_array + idxmm - m;
          al = &alm[2 * (idxmm + 2)];            /* a(m+2,m) */
          while (l < lmax)
            {
              plm    = (al[0] * x) * plmp1 + al[1] * plm;
              dplm   = al[0] * (x * dplmp1 - usq * plmp1) + al[1] * dplm;

              plmp1  = (al[2] * x) * plm + al[3] * plmp1;
              dplmp1 = al[2] * (x * dplm - usq * plm) + al[3] * dplmp1;

              ptr[l] = plm; ptr[l+1] = plmp1;
              dptr[l] = dplm; dptr[l+1] = dplmp1;

              l += 2; al += 4;
            }

          if (l == lmax)
            {
              ptr[l] = (al[0] * x) * plmp1 + al[1] * plm;
              dptr[l] = al[0] * (x * dplmp1 - usq * plmp1) + al[1] * dplm;
            }

          if (m < mmax)
            {
              /* compute P(m+1,m+1) = d_{m+1} * u * P(m,m) */
              idxmm += Lp1 - m; /* idx(m+1,m+1) = idx(m,m) + L + 1 - m */
              dpmm = dl[m+1] * u * (x * pmm + dpmm);
              pmm *= dl[m+1] * u;
              result_array[idxmm] = pmm;
              result_deriv_array[idxmm] = dpmm;
            }
        }

      return GSL_SUCCESS;
    }
}

#if 0
/*
gsl_sf_alf_theta_deriv2_array()
  Compute array of associated Legendre functions and their
first and second theta derivatives at a given point x

Inputs: norm                - ALF normalization
        lmax                - maximum degree
        mmax                - maximum order
        x                   - input point
        result_array        - (input/output) output array of ALFs
                              Plm(x)            for m = 0
                              Plm(x)/sin(theta) for m > 0
        result_deriv_array  - (output) output array of ALF derivatives
                              d/dtheta Plm(x)
        result_deriv2_array - (output) output array of second ALF derivatives
                              d^2/dtheta^2 Plm(x)

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) result_array, result_deriv_array, and result_deriv2_array are
indexed by gsl_sf_alf_array_index(l,m,lmax),

result_array[index(l,0,lmax)]        = Pl0(x)
result_array[index(l,m,lmax)]        = Plm(x)/sin(theta), m > 0
result_deriv_array[index(l,m,lmax)]  = dPlm(x)/dtheta
result_deriv2_array[index(l,m,lmax)] = d^2 Plm(x)/dtheta^2
*/

int
gsl_sf_alf_theta_deriv2_array(const gsl_sf_alf_t norm, const size_t lmax, const size_t mmax,
                              const double x, double result_array[], double result_deriv_array[],
                              double result_deriv2_array[])
{
  int status;

  status = gsl_sf_alf_theta_deriv_array(lmax, mmax, x, result_array, result_deriv_array);
  if (status)
    return status;

  status = gsl_sf_alf_theta_derivk_array(norm, lmax, mmax, result_array,
                                         result_deriv_array, result_deriv2_array);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
gsl_sf_alf_theta_derivk_array()
  Compute the kth derivative of Plm(x) with respect to theta
using Eq. 14 of [2]

Inputs: norm         - ALF normalization
        lmax         - maximum degree
        mmax         - maximum order
        Plm          - array initialized by gsl_sf_alf_precompute()
        input_array  - (input) array of d^{k-1}/dtheta^{k-1} Plm(x), length nlm
        output_array - (output) d^k/dtheta^k Plm(x), length nlm

Notes:
1) Plm can equal input_array to calculate first derivatives
*/

int
gsl_sf_alf_theta_derivk_array(const gsl_sf_alf_t norm,
                              const size_t lmax,
                              const size_t mmax,
                              const double Plm[],
                              const double input_array[],
                              double output_array[])
{
  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else
    {
      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      const double * alm = &Plm[nlm];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      const double * sqrts = dl + lmax + 1;
      const double csfac = dl[0];
      const size_t Lp1 = lmax + 1;
      size_t l, m;
      size_t idxm, idxmm1, idxmp1;

      /* d^k/dtheta^k P(0,0) = 0 */
      output_array[0] = 0.0;

      if (lmax == 0)
        return GSL_SUCCESS;

      if (norm == GSL_SF_ALF_NONE)
        {
          idxmm1 = 0;      /* idx(0,0) */
          idxm = lmax + 1; /* idx(1,1) */
          for (l = 1; l <= lmax; ++l)
            {
              /* Eq. 13a: compute d^k/dtheta^k P(l,0) = -d^{k-1}/dtheta^{k-1} P(l,1) */
              output_array[l] = -csfac * input_array[lmax + l];

              /* Eq. 13b: compute d^k/dtheta^k P(l,l) = l d^{k-1}/dtheta^{k-1} P(l,l-1) */
              output_array[idxm] = csfac * l * input_array[idxmm1 + 1];

              idxmm1 = idxm;
              idxm += Lp1 - l; /* idx(l+1,l+1) */
            }

          idxmm1 = 0;            /* idx(0,0) */
          idxm = lmax + 1;       /* idx(1,1) */
          idxmp1 = 2 * lmax + 1; /* idx(2,2) */
          for (m = 1; m < mmax; ++m)
            {
              size_t offset = 1;

              for (l = m + 1; l <= lmax; ++l)
                {
                  output_array[idxm + offset] = 0.5 * csfac *
                                                ((l + m) * (l - m + 1.0) * input_array[idxmm1 + offset + 1] -
                                                 input_array[idxmp1 + offset - 1]);
                  ++offset;
                }

              idxmm1 = idxm;
              idxm = idxmp1;
              idxmp1 += Lp1 - m - 1; /* idx(m+1,m+1) */
            }
        }
      else if ((norm == GSL_SF_ALF_SCHMIDT) ||
               (norm == GSL_SF_ALF_FOURPI))
        {
          /* d^k/dtheta^k P(1,0) = -d^{k-1}/dtheta^{k-1} P(1,1) */
          output_array[1] = -csfac * input_array[lmax + 1];

          /* d^k/dtheta^k P(1,1) = d^{k-1}/dtheta^{k-1} P(1,0) */
          output_array[lmax + 1] = csfac * input_array[1];

          /* compute m=0 terms */
          idxmm1 = lmax + 1;   /* idx(1,1) */
          idxm = 2 * lmax + 1; /* idx(2,2) */
          for (l = 2; l <= lmax; ++l)
            {
              /* d^k/dtheta^k P(l,0) = -sqrt(l(l+1)) d^{k-1}/dtheta^{k-1} P(l,1) */
              output_array[l] = -csfac * (sqrts[l] / M_SQRT2) * sqrts[l + 1] * input_array[lmax + l];

              /* d^k/dtheta^k P(l,l) = sqrt(l/2) d^{k-1}/dtheta^{k-1} P(l,l-1) */
              output_array[idxm] = csfac * (sqrts[l] / M_SQRT2) * input_array[idxmm1 + 1];

              idxmm1 = idxm;
              idxm += Lp1 - l; /* idx(l+1,l+1) */
            }

          /* compute m=1 terms */
          idxmm1 = 0;            /* idx(0,0) */
          idxm = lmax + 1;       /* idx(1,1) */
          idxmp1 = 2 * lmax + 1; /* idx(2,2) */
          for (l = 2; l <= lmax; ++l)
            {
              output_array[idxm + l - 1] = csfac *
                                          ((sqrts[l] / M_SQRT2) * sqrts[l + 1] * input_array[idxmm1 + l] -
                                           0.5 * sqrts[l + 2] * sqrts[l - 1] * input_array[idxmp1 + l - 2]);
            }

          /* compute m>1 terms */
          idxmm1 = lmax + 1;   /* idx(1,1) */
          idxm = 2 * lmax + 1; /* idx(2,2) */
          idxmp1 = 3 * lmax;   /* idx(3,3) */
          for (m = 2; m < mmax; ++m)
            {
              size_t offset = 1;

              for (l = m + 1; l <= lmax; ++l)
                {
                  output_array[idxm + offset] = 0.5 * csfac *
                                                (sqrts[l + m] * sqrts[l - m + 1] * input_array[idxmm1 + offset + 1] -
                                                 sqrts[l + m + 1] * sqrts[l - m] * input_array[idxmp1 + offset - 1]);
                  ++offset;
                }

              idxmm1 = idxm;
              idxm = idxmp1;
              idxmp1 += Lp1 - m - 1; /* idx(m+1,m+1) */
            }
        }
      else
        {
          idxmm1 = 0;      /* idx(0,0) */
          idxm = lmax + 1; /* idx(1,1) */
          for (l = 1; l <= lmax; ++l)
            {
              /* d^k/dtheta^k P(l,0) = -sqrt(l(l+1)) d^{k-1}/dtheta^{k-1} P(l,1) */
              output_array[l] = -csfac * sqrts[l] * sqrts[l + 1] * input_array[lmax + l];

              /* d^k/dtheta^k P(l,l) = sqrt(l/2) d^{k-1}/dtheta^{k-1} P(l,l-1) */
              output_array[idxm] = csfac * (sqrts[l] / M_SQRT2) * input_array[idxmm1 + 1];

              idxmm1 = idxm;
              idxm += Lp1 - l; /* idx(l+1,l+1) */
            }

          idxmm1 = 0;            /* idx(0,0) */
          idxm = lmax + 1;       /* idx(1,1) */
          idxmp1 = 2 * lmax + 1; /* idx(2,2) */
          for (m = 1; m < mmax; ++m)
            {
              size_t offset = 1;

              for (l = m + 1; l <= lmax; ++l)
                {
                  output_array[idxm + offset] = 0.5 * csfac *
                                                (sqrts[l + m] * sqrts[l - m + 1] * input_array[idxmm1 + offset + 1] -
                                                 sqrts[l + m + 1] * sqrts[l - m] * input_array[idxmp1 + offset - 1]);
                  ++offset;
                }

              idxmm1 = idxm;
              idxm = idxmp1;
              idxmp1 += Lp1 - m - 1; /* idx(m+1,m+1) */
            }
        }

      return GSL_SUCCESS;
    }
}
#endif

#if 0

/*
gsl_sf_alf_deriv2_alt_arrayx()
  Compute array of associated Legendre functions and their
first and second derivatives at a given point x

Inputs: norm                - ALF normalization
        lmax                - maximum degree
        x                   - input point, cos(theta)
        result_array        - (input/output) output array of ALFs Plm(x)
        result_deriv_array  - (output) d/dx Plm(x), length nlm
        result_deriv2_array - (output) d^2/dx^2 Plm(x), length nlm

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) result_deriv_array and result_deriv2_array must be of length at least nlm

3) output arrays are indexed by gsl_sf_alf_array_index(l,m)
*/

int
gsl_sf_alf_deriv2_arrayx(const gsl_sf_alf_t norm,
                              const size_t lmax,
                              const double x,
                              double result_array[],
                              double result_deriv_array[],
                              double result_deriv2_array[])
{
  if (x == -1.0 || x == 1.0)
    {
      GSL_ERROR ("x cannot equal -1 or +1 for derivative computation", GSL_EDOM);
    }
  else
    {
      int status;
      const size_t nlm = gsl_sf_alf_nlm(lmax);
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double uinv = 1.0 / u;
      const double uinv2 = uinv * uinv;
      size_t i;

      /* compute Plm, d/dtheta Plm and d^2/dtheta^2 Plm */
      status = gsl_sf_alf_deriv2_alt_arrayx(norm, lmax, x,
                                                 result_array,
                                                 result_deriv_array,
                                                 result_deriv2_array);
      if (status)
        return status;

      for (i = 0; i < nlm; ++i)
        {
          /* d/dx Plm = -1/sin(theta) d/dtheta Plm */
          result_deriv_array[i] *= -uinv;

          /* d^2/dx^2 Plm = 1/sin^2(theta) [ d^2/dtheta^2 Plm + x d/dx Plm ] */
          result_deriv2_array[i] = uinv2 * (result_deriv2_array[i] + x * result_deriv_array[i]);
        }

      return GSL_SUCCESS;
    }
}

#endif
