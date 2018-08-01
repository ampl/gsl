/* statistics/Qn_source.c
 * 
 * Copyright (C) 2018 Patrick Alken
 * Copyright (C) 2005, 2006, 2007 Martin Maechler, ETH Zurich
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

/* This is a merge of the C version of original files  qn.f and sn.f,
 * translated by f2c (version 20010821).	       ====	====
 * and then by f2c-clean,v 1.9 2000/01/13 13:46:53
 * and further clean-edited manually by Martin Maechler.
 *
 * Further added interface functions to be called via .C() from R or S-plus
 * Note that Peter Rousseeuw has explicitely given permission to
 * use his code under the GPL for the R project.
 */

/* Original comments by the authors of the Fortran original code,
 * (merged for Qn & Sn in one file by M.M.):

   This file contains fortran functions for two new robust estimators
   of scale denoted as Qn and Sn, decribed in Rousseeuw and Croux (1993).
   These estimators have a high breakdown point and a bounded influence
   function. The implementation given here is very fast (running in
   O(n logn) time) and needs little storage space.

	Rousseeuw, P.J. and Croux, C. (1993)
	Alternatives to the Median Absolute Deviation",
	Journal of the American Statistical Association, Vol. 88, 1273-1283.

   For both estimators, implementations in the pascal language can be
   obtained from the original authors.

   This software may be used and copied freely for scientific
   and/or non-commercial purposes, provided reference is made
   to the abovementioned paper.

Note by MM: We have explicit permission from P.Rousseeuw to
licence it under the GNU Public Licence.
*/

#ifndef int64_t
#define int64_t long int
#endif

static BASE FUNCTION(Qn,whimed)(BASE * a, int * w, int n, BASE * a_cand, BASE * a_srt, int * w_cand);

/*
gsl_stats_Qn0_from_sorted_data()
  Efficient algorithm for the scale estimator:

    Q_n0 = { |x_i - x_j|; i<j }_(k) [ = Qn without scaling ]

i.e. the k-th order statistic of the |x_i - x_j|, where:

k = (floor(n/2) + 1 choose 2)

Inputs: sorted_data - sorted array containing the observations
        stride      - stride
        n           - length of 'sorted_data'
        work        - workspace of length 3n of type BASE
        work_int    - workspace of length 5n of type int

Return: Q_n statistic (without scale/correction factor); same type as input data
*/

BASE
FUNCTION(gsl_stats,Qn0_from_sorted_data) (const BASE sorted_data[],
                                          const size_t stride,
                                          const size_t n,
                                          BASE work[],
                                          int work_int[])
{
  const int ni = (int) n;
  BASE * a_srt = &work[n];
  BASE * a_cand = &work[2*n];

  int *left = &work_int[0];
  int *right = &work_int[n];
  int *p = &work_int[2*n];
  int *q = &work_int[3*n];
  int *weight = &work_int[4*n];

  BASE trial = (BASE) 0.0;
  int found = 0;

  int h, i, j, jh;

  /* following should be `long long int' : they can be of order n^2 */
  int64_t k, knew, nl,nr, sump, sumq;

  /* check for quick return */
  if (n < 2)
    return ((BASE) 0.0);

  h = n / 2 + 1;
  k = (int64_t)h * (h - 1) / 2;

  for (i = 0; i < ni; ++i)
    {
      left[i] = ni - i + 1;
      right[i] = (i <= h) ? ni : ni - (i - h);

      /* the n - (i-h) is from the paper; original code had `n' */
    }

  nl = (int64_t)n * (n + 1) / 2;
  nr = (int64_t)n * n;
  knew = k + nl;/* = k + (n+1 \over 2) */

/* L200: */
  while (!found && nr - nl > ni)
    {
      j = 0;
      /* Truncation to float : try to make sure that the same values are got later (guard bits !) */
      for (i = 1; i < ni; ++i)
        {
          if (left[i] <= right[i])
            {
              weight[j] = right[i] - left[i] + 1;
              jh = left[i] + weight[j] / 2;
              work[j] = sorted_data[i * stride] - sorted_data[(ni - jh) * stride];
              ++j;
            }
        }

      trial = FUNCTION(Qn,whimed)(work, weight, j, a_cand, a_srt, /*iw_cand*/ p);

      j = 0;
      for (i = ni - 1; i >= 0; --i)
        {
          while (j < ni && ((double)(sorted_data[i * stride] - sorted_data[(ni - j - 1) * stride])) < trial)
            ++j;

          p[i] = j;
        }

      j = ni + 1;
      for (i = 0; i < ni; ++i)
        {
          while ((double)(sorted_data[i * stride] - sorted_data[(ni - j + 1) * stride]) > trial)
            --j;

          q[i] = j;
        }

      sump = 0;
      sumq = 0;

      for (i = 0; i < ni; ++i)
        {
          sump += p[i];
          sumq += q[i] - 1;
        }

      if (knew <= sump)
        {
          for (i = 0; i < ni; ++i)
            right[i] = p[i];

          nr = sump;
        }
      else if (knew > sumq)
        {
          for (i = 0; i < ni; ++i)
            left[i] = q[i];

          nl = sumq;
        }
      else /* sump < knew <= sumq */
        {
          found = 1;
        }
    } /* while */

  if (found)
    {
      return trial;
    }
  else
    {
      j = 0;
      for (i = 1; i < ni; ++i)
        {
          int jj;

          for (jj = left[i]; jj <= right[i]; ++jj)
            {
              work[j] = sorted_data[i * stride] - sorted_data[(ni - jj) * stride];
              j++;
            }/* j will be = sum_{i=2}^n (right[i] - left[i] + 1)_{+}  */
        }

      /* return pull(work, j - 1, knew - nl)	: */
      knew -= (nl + 1); /* -1: 0-indexing */

      /* sort work array */
      TYPE (gsl_sort) (work, 1, j);

      return (work[knew]);
    }
}

/*
gsl_stats_Qn_from_sorted_data()
  Efficient algorithm for the scale estimator:

    Q_n = 2.219 * d_n * { |x_i - x_j|; i<j }_(k)

with:

k = (floor(n/2) + 1 choose 2)

and d_n is a correction factor for finite sample bias

Inputs: sorted_data - sorted array containing the observations
        stride      - stride
        n           - length of 'sorted_data'
        work        - workspace of length 3n of type BASE
        work_int    - workspace of length 5n of type int

Return: Q_n statistic
*/

double
FUNCTION(gsl_stats,Qn_from_sorted_data) (const BASE sorted_data[],
                                         const size_t stride,
                                         const size_t n,
                                         BASE work[],
                                         int work_int[])
{
  const double scale = 2.21914; /* asymptotic consistency for sigma^2 */
  double Qn0 = (double) FUNCTION(gsl_stats,Qn0_from_sorted_data)(sorted_data, stride, n, work, work_int);
  double dn = 1.0;
  double Qn;

  /* this correction factor deviates from the original paper Croux and Rousseeuw, 1992, and
   * comes from the 'robustbase' R package */
  if (n <= 12)
    {
      if (n == 2) dn = .399356;
      else if (n == 3) dn = .99365;
      else if (n == 4) dn = .51321;
      else if (n == 5) dn = .84401;
      else if (n == 6) dn = .61220;
      else if (n == 7) dn = .85877;
      else if (n == 8) dn = .66993;
      else if (n == 9) dn = .87344;
      else if (n == 10) dn = .72014;
      else if (n == 11) dn = .88906;
      else if (n == 12) dn = .75743;
    }
  else
    {
      if (n % 2 == 1) /* n odd */
        dn = 1.60188 + (-2.1284 - 5.172 / n) / n;
      else            /* n even */
        dn = 3.67561 + (1.9654 + (6.987 - 77.0 / n) / n) / n;

      dn = 1.0 / (dn / (double)n + 1.0);
    }

  Qn = scale * dn * Qn0;

  return Qn;
}

/*
  Algorithm to compute the weighted high median in O(n) time.

  The whimed is defined as the smallest a[j] such that the sum
  of the weights of all a[i] <= a[j] is strictly greater than
  half of the total weight.

  Arguments:

  a: double array containing the observations
  n: number of observations
  w: array of (int/double) weights of the observations.
*/

static BASE
FUNCTION(Qn,whimed)(BASE * a, int * w, int n, BASE * a_cand, BASE * a_srt, int * w_cand)
{
  int n2, i, kcand;
  /* sum of weights: `int' do overflow when  n ~>= 1e5 */
  int64_t wleft, wmid, wright, w_tot, wrest;

  BASE trial;

  w_tot = 0;
  for (i = 0; i < n; ++i)
    w_tot += w[i];

  wrest = 0;

  /* REPEAT : */
  do
    {
      for (i = 0; i < n; ++i)
        a_srt[i] = a[i];

      n2 = n/2; /* =^= n/2 +1 with 0-indexing */
#if 0
      rPsort(a_srt, n, n2);
#else
      TYPE (gsl_sort) (a_srt, 1, n);
#endif
      trial = a_srt[n2];

      wleft = 0;
      wmid = 0;
      wright = 0;

      for (i = 0; i < n; ++i)
        {
          if (a[i] < trial)
            wleft += w[i];
          else if (a[i] > trial)
            wright += w[i];
          else
            wmid += w[i];
        }

      /* wleft = sum_{i; a[i]	 < trial}  w[i]
       * wmid	 = sum_{i; a[i] == trial}  w[i] at least one 'i' since trial is one a[]!
       * wright= sum_{i; a[i]	 > trial}  w[i]
       */
      kcand = 0;
      if (2 * (wrest + wleft) > w_tot)
        {
          for (i = 0; i < n; ++i)
            {
              if (a[i] < trial)
                {
                  a_cand[kcand] = a[i];
                  w_cand[kcand] = w[i];
                  ++kcand;
                }
            }
        }
      else if (2 * (wrest + wleft + wmid) <= w_tot)
        {
          for (i = 0; i < n; ++i)
            {
              if (a[i] > trial)
                {
                  a_cand[kcand] = a[i];
                  w_cand[kcand] = w[i];
                  ++kcand;
                }
            }

          wrest += wleft + wmid;
        }
      else
        {
          return trial;
        }

      n = kcand;
      for (i = 0; i < n; ++i)
        {
          a[i] = a_cand[i];
          w[i] = w_cand[i];
        }
    } while(1);
}
