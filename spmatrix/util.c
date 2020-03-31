/* spmatrix/util.c
 * 
 * Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Patrick Alken
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>

/*
gsl_spmatrix_cumsum()

Compute the cumulative sum:

p[j] = Sum_{k=0...j-1} c[k]

0 <= j < n + 1

Alternatively,
p[0] = 0
p[j] = p[j - 1] + c[j - 1]

Inputs: n - length of input array
        c - (input/output) array of size n + 1
            on input, contains the n values c[k]
            on output, contains the n + 1 values p[j]

Return: success or error
*/

void
gsl_spmatrix_cumsum(const size_t n, int * c)
{
  int sum = 0;
  size_t k;

  for (k = 0; k < n; ++k)
    {
      int ck = c[k];
      c[k] = sum;
      sum += ck;
    }

  c[n] = sum;
}
