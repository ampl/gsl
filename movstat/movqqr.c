/* movstat/movqqr.c
 *
 * Compute moving q-quantile range
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
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

/*
gsl_movstat_qqr()
  Apply a moving q-quantile range to an input vector

Inputs: endtype - how to handle end points
        x       - input vector, size n
        q       - quantile \in [0,0.5]
        xqqr    - (output) vector of q-quantile ranges of x, size n
                  xqqr_i = q-quantile range of i-th window:
        w       - workspace
*/

int
gsl_movstat_qqr(const gsl_movstat_end_t endtype, const gsl_vector * x, const double q,
                gsl_vector * xqqr, gsl_movstat_workspace * w)
{
  if (x->size != xqqr->size)
    {
      GSL_ERROR("x and xqqr vectors must have same length", GSL_EBADLEN);
    }
  else if (q < 0.0 || q > 0.5)
    {
      GSL_ERROR("q must be between 0 and 0.5", GSL_EDOM);
    }
  else
    {
      double qq = q;
      int status = gsl_movstat_apply_accum(endtype, x, gsl_movstat_accum_qqr, (void *) &qq, xqqr, NULL, w);
      return status;
    }
}
