/* movstat/movmad.c
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
gsl_movstat_mad()
  Apply a moving MAD to an input vector (with scale factor to make an
unbiased estimate of sigma)

Inputs: endtype - how to handle end points
        x       - input vector, size n
        xmedian - (output) vector of median values of x, size n
                  xmedian_i = median of window centered on x_i
        xmad    - (output) vector of estimated standard deviations of x, size n
                  xmad_i = MAD of i-th window: 1.4826 * median(|x_i - xmedian_i|)
        w       - workspace
*/

int
gsl_movstat_mad(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * xmedian, gsl_vector * xmad,
                gsl_movstat_workspace * w)
{
  if (x->size != xmedian->size)
    {
      GSL_ERROR("x and xmedian vectors must have same length", GSL_EBADLEN);
    }
  else if (x->size != xmad->size)
    {
      GSL_ERROR("x and xmad vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      double scale = 1.482602218505602;
      int status = gsl_movstat_apply_accum(endtype, x, gsl_movstat_accum_mad, &scale, xmedian, xmad, w);
      return status;
    }
}

/*
gsl_movstat_mad0()
  Apply a moving MAD to an input vector (without scale factor to make an
unbiased estimate of sigma)

Inputs: endtype - how to handle end points
        x       - input vector, size n
        xmedian - (output) vector of median values of x, size n
                  xmedian_i = median of window centered on x_i
        xmad    - (output) vector of estimated standard deviations of x, size n
                  xmad_i = MAD of i-th window: median(|x_i - xmedian_i|)
        w       - workspace
*/

int
gsl_movstat_mad0(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * xmedian, gsl_vector * xmad,
                 gsl_movstat_workspace * w)
{
  if (x->size != xmedian->size)
    {
      GSL_ERROR("x and xmedian vectors must have same length", GSL_EBADLEN);
    }
  else if (x->size != xmad->size)
    {
      GSL_ERROR("x and xmad vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      double scale = 1.0;
      int status = gsl_movstat_apply_accum(endtype, x, gsl_movstat_accum_mad, &scale, xmedian, xmad, w);
      return status;
    }
}
