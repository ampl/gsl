/* statistics/gastwirth_source.c
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
FUNCTION(gsl_stats,gastwirth_from_sorted_data) (const BASE sorted_data[],
                                                const size_t stride,
                                                const size_t n)
{
  if (n == 0)
    {
      return 0.0;
    }
  else
    {
      double a = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted_data, stride, n, 1.0 / 3.0);
      double b = FUNCTION(gsl_stats,median_from_sorted_data)(sorted_data, stride, n);
      double c = FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted_data, stride, n, 2.0 / 3.0);
      double gastwirth = 0.3 * a + 0.4 * b + 0.3 * c;

      return gastwirth;
    }
}
