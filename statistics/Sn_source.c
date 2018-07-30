/* statistics/Sn_source.c
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

/*
gsl_stats_Sn0_from_sorted_data()
  Efficient algorithm for the scale estimator:

    S_n0 = LOMED_{i} HIMED_{i} |x_i - x_j|

which can equivalently be written as

    S_n0 = LOMED_{i} LOMED_{j != i} |x_i - x_j|

Inputs: sorted_data - sorted array containing the observations
        stride      - stride
        n           - length of 'sorted_data'
        work        - workspace of length n
                      work[i] := LOMED_{j != i} | x_i - x_j |

Return: S_n statistic (without scale/correction factor)
*/

BASE
FUNCTION(gsl_stats,Sn0_from_sorted_data) (const BASE sorted_data[],
                                          const size_t stride,
                                          const size_t n,
                                          BASE work[])
{
  /* Local variables */
  double medA, medB;
  int i, diff, half, Amin, Amax, even, length;
  int leftA, leftB, nA, nB, tryA, tryB, rightA, rightB;
  int np1_2 = (n + 1) / 2;

  work[0] = sorted_data[n / 2 * stride] - sorted_data[0];

  /* first half for() loop : */
  for (i = 2; i <= np1_2; ++i)
    {
      nA = i - 1;
      nB = n - i;
      diff = nB - nA;
      leftA = leftB	= 1;
      rightA = rightB = nB;
      Amin = diff / 2 + 1;
      Amax = diff / 2 + nA;

      while (leftA < rightA)
        {
          length = rightA - leftA + 1;
          even = 1 - length % 2;
          half = (length - 1) / 2;
          tryA = leftA + half;
          tryB = leftB + half;
          if (tryA < Amin)
            {
              rightB = tryB;
              leftA = tryA + even;
            }
          else
            {
              if (tryA > Amax)
                {
                  rightA = tryA;
                  leftB = tryB + even;
                }
              else
                {
                  medA = sorted_data[(i - 1) * stride] - sorted_data[(i - tryA + Amin - 2) * stride];
                  medB = sorted_data[(tryB + i - 1) * stride] - sorted_data[(i - 1) * stride];
                  if (medA >= medB)
                    {
                      rightA = tryA;
                      leftB = tryB + even;
                    }
                  else
                    {
                      rightB = tryB;
                      leftA = tryA + even;
                    }
                }
            }
        } /* while */

      if (leftA > Amax)
        {
          work[i - 1] = sorted_data[(leftB + i - 1) * stride] - sorted_data[(i - 1) * stride];
        }
      else
        {
          medA = sorted_data[(i - 1) * stride] - sorted_data[(i - leftA + Amin - 2) * stride];
          medB = sorted_data[(leftB + i - 1) * stride] - sorted_data[(i - 1) * stride];
          work[i - 1] = GSL_MIN(medA, medB);
        }
    }

  /* second half for() loop : */
  for (i = np1_2 + 1; i <= (int) n - 1; ++i)
    {
      nA = n - i;
      nB = i - 1;
      diff = nB - nA;
      leftA  = leftB	= 1;
      rightA = rightB = nB;
      Amin = diff / 2 + 1;
      Amax = diff / 2 + nA;

      while (leftA < rightA)
        {
          length = rightA - leftA + 1;
          even = 1 - length % 2;
          half = (length - 1) / 2;
          tryA = leftA + half;
          tryB = leftB + half;

          if (tryA < Amin)
            {
              rightB = tryB;
              leftA = tryA + even;
            }
          else
            {
              if (tryA > Amax)
                {
                  rightA = tryA;
                  leftB = tryB + even;
                }
              else
                {
                  medA = sorted_data[(i + tryA - Amin) * stride] - sorted_data[(i - 1) * stride];
                  medB = sorted_data[(i - 1) * stride] - sorted_data[(i - tryB - 1) * stride];
                  if (medA >= medB)
                    {
                      rightA = tryA;
                      leftB = tryB + even;
                    }
                  else
                    {
                      rightB = tryB;
                      leftA = tryA + even;
                    }
                }
            }
        } /* while */

      if (leftA > Amax)
        {
          work[i - 1] = sorted_data[(i - 1) * stride] - sorted_data[(i - leftB - 1) * stride];
        }
      else
        {
          medA = sorted_data[(i + leftA - Amin) * stride] - sorted_data[(i - 1) * stride];
          medB = sorted_data[(i - 1) * stride] - sorted_data[(i - leftB - 1) * stride];
          work[i - 1] = GSL_MIN(medA, medB);
        }
    }

  work[n - 1] = sorted_data[(n - 1) * stride] - sorted_data[(np1_2 - 1) * stride];

  /* sort work array */
  TYPE (gsl_sort) (work, 1, n);

  return work[np1_2 - 1];
}

/*
gsl_stats_Sn_from_sorted_data()
  Efficient algorithm for the scale estimator:

    S_n = 1.1926 * c_n LOMED_{i} HIMED_{i} |x_i - x_j|

which can equivalently be written as

    S_n = 1.1926 * c_n * LOMED_{i} LOMED_{j != i} |x_i - x_j|

and c_n is a correction factor for small sample bias

Inputs: sorted_data - sorted array containing the observations
        stride      - stride
        n           - length of 'sorted_data'
        work        - workspace of length n
                      work[i] := LOMED_{j != i} | x_i - x_j |

Return: S_n statistic
*/

double
FUNCTION(gsl_stats,Sn_from_sorted_data) (const BASE sorted_data[],
                                         const size_t stride,
                                         const size_t n,
                                         BASE work[])
{
  const double scale = 1.1926; /* asymptotic consistency for sigma^2 */
  double Sn0 = (double) FUNCTION(gsl_stats,Sn0_from_sorted_data)(sorted_data, stride, n, work);
  double cn = 1.0;
  double Sn;

  /* determine correction factor for finite sample bias */
  if (n <= 9)
    {
      if (n == 2) cn = 0.743;
      else if (n == 3) cn = 1.851;
      else if (n == 4) cn = 0.954;
      else if (n == 5) cn = 1.351;
      else if (n == 6) cn = 0.993;
      else if (n == 7) cn = 1.198;
      else if (n == 8) cn = 1.005;
      else if (n == 9) cn = 1.131;
    }
  else if (n % 2 == 1) /* n odd, >= 11 */
    {
      cn = (double) n / (n - 0.9);
    }

  Sn = scale * cn * Sn0;

  return Sn;
}
