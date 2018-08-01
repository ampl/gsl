/* movstat/movSn.c
 *
 * Compute moving "S_n" statistic from Croux and Rousseeuw, 1992
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
gsl_movstat_Sn()
  Calculate moving S_n statistic for input vector

Inputs: endtype - how to handle end points
        x       - input vector, size n
        xscale  - (output) vector of "S_n" statistics, size n
                  xscale_i = (S_n)_i for i-th window:
        w       - workspace
*/

int
gsl_movstat_Sn(const gsl_movstat_end_t endtype, const gsl_vector * x,
               gsl_vector * xscale, gsl_movstat_workspace * w)
{
  int status = gsl_movstat_apply_accum(endtype, x, gsl_movstat_accum_Sn, NULL, xscale, NULL, w);
  return status;
}
