/* movstat/movminmax.c
 *
 * Routines related to a moving window min/max
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_movstat.h>

/*
gsl_movstat_minmax()
  Apply minmax filter to input vector

Inputs: endtype - end point handling criteria
        x       - input vector, size n
        y_min   - output vector of minimum values, size n
        y_max   - output vector of maximum values, size n
        w       - workspace
*/

int
gsl_movstat_minmax(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y_min, gsl_vector * y_max, gsl_movstat_workspace * w)
{
  int status = gsl_movstat_apply_accum(endtype, x, gsl_movstat_accum_minmax, NULL, y_min, y_max, w);
  return status;
}

/*
gsl_movstat_min()
  Apply minimum filter to input vector

Inputs: endtype - end point handling criteria
        x       - input vector, size n
        y       - output vector of minimum values, size n
        w       - workspace
*/

int
gsl_movstat_min(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)
{
  int status = gsl_movstat_apply_accum(endtype, x, gsl_movstat_accum_min, NULL, y, NULL, w);
  return status;
}

/*
gsl_movstat_max()
  Apply maximum filter to input vector

Inputs: endtype - end point handling criteria
        x       - input vector, size n
        y       - output vector of maximum values, size n
        w       - workspace
*/

int
gsl_movstat_max(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)
{
  int status = gsl_movstat_apply_accum(endtype, x, gsl_movstat_accum_max, NULL, y, NULL, w);
  return status;
}
