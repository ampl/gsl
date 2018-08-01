/* statistics/mad_source.c
 * 
 * Copyright (C) 2018 Patrick Alken
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

/*
gsl_stats_mad0()
  Compute median absolute deviation

Inputs: data    - array containing the observations
        stride  - stride
        n       - length of 'data'
        work    - workspace of length n

Return: MAD statistic (without scale/correction factor)
*/

double
FUNCTION(gsl_stats,mad0) (const BASE data[],
                          const size_t stride,
                          const size_t n,
                          double work[])
{
  double median, mad;
  size_t i;

  /* copy input data to work */
  for (i = 0; i < n; ++i)
    work[i] = (double) data[i * stride];

  /* compute median of input data using double version */
  median = gsl_stats_median(work, 1, n);

  /* compute absolute deviations from median */
  for (i = 0; i < n; ++i)
    work[i] = fabs((double) data[i * stride] - median);

  mad = gsl_stats_median(work, 1, n);

  return mad;
}

double
FUNCTION(gsl_stats,mad) (const BASE data[],
                         const size_t stride,
                         const size_t n,
                         double work[])
{
  double mad0 = FUNCTION(gsl_stats,mad0)(data, stride, n, work);
  double mad = 1.482602218505602 * mad0;
  return mad;
}
