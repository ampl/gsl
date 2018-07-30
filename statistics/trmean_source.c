/* statistics/trmean_source.c
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

double
FUNCTION (gsl_stats, trmean_from_sorted_data) (const double trim, const BASE sorted_data[],
                                               const size_t stride, const size_t size)
{
  if (trim >= 0.5)
    {
      return FUNCTION(gsl_stats, median_from_sorted_data)(sorted_data, stride, size);
    }
  else
    {
      size_t ilow = (size_t) floor(trim * size);
      size_t ihigh = size - ilow - 1;
      double mean = 0.0;
      double k = 0.0;
      size_t i;

      /* compute mean of middle samples in [ilow,ihigh] */
      for (i = ilow; i <= ihigh; ++i)
        {
          double delta = sorted_data[i * stride] - mean;
          k += 1.0;
          mean += delta / k;
        }

      return mean;
    }
}
