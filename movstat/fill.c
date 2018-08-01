/* movstat/fill.c
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
 
#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>

/*
gsl_movstat_fill()
  Fill window for sample 'idx' from x using given end conditions

Inputs: endtype     - how to handle end points
        x           - input vector, length n
        idx         - index of center sample in window in [0,n-1]
        H           - number of samples left of center to include
        J           - number of samples right of center to include
        window      - (output) window of samples centered on x_{idx},
                      W_{idx}^{H,J}, length H + J + 1

Return: size of filled window (<= H + J + 1)
*/

size_t
gsl_movstat_fill(const gsl_movstat_end_t endtype, const gsl_vector * x, const size_t idx,
                 const size_t H, const size_t J, double * window)
{
  if (idx >= x->size)
    {
      GSL_ERROR_VAL ("window center index must be between 0 and n - 1", GSL_EDOM, 0);
    }
  else
    {
      const int n = (int) x->size;
      const int iidx = (int) idx;
      const int iH = (int) H;
      const int iJ = (int) J;
      int idx1, idx2, j;
      size_t window_size;

      if (endtype == GSL_MOVSTAT_END_TRUNCATE)
        {
          idx1 = GSL_MAX(iidx - iH, 0);
          idx2 = GSL_MIN(iidx + iJ, n - 1);
        }
      else
        {
          idx1 = iidx - iH;
          idx2 = iidx + iJ;
        }

      window_size = (size_t) (idx2 - idx1 + 1);

      /* fill sliding window */
      for (j = idx1; j <= idx2; ++j)
        {
          int widx = j - idx1;

          if (j < 0)
            {
              /* initial condition */
              if (endtype == GSL_MOVSTAT_END_PADZERO)
                window[widx] = 0.0;
              else if (endtype == GSL_MOVSTAT_END_PADVALUE)
                window[widx] = gsl_vector_get(x, 0);
            }
          else if (j >= n)
            {
              if (endtype == GSL_MOVSTAT_END_PADZERO)
                window[widx] = 0.0;
              else if (endtype == GSL_MOVSTAT_END_PADVALUE)
                window[widx] = gsl_vector_get(x, n - 1);
            }
          else
            {
              window[widx] = gsl_vector_get(x, j);
            }
        }

      return window_size;
    }
}
