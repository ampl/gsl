/* filter/median.c
 *
 * Contains routines related to the standard median filter
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_errno.h>

gsl_filter_median_workspace *
gsl_filter_median_alloc(const size_t K)
{
  gsl_filter_median_workspace *w;
  size_t H = K / 2;

  w = calloc(1, sizeof(gsl_filter_median_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->movstat_workspace_p = gsl_movstat_alloc(2*H + 1);
  if (w->movstat_workspace_p == NULL)
    {
      gsl_filter_median_free(w);
      GSL_ERROR_NULL ("failed to allocate space for movstat workspace", GSL_ENOMEM);
    }

  return w;
}

void
gsl_filter_median_free(gsl_filter_median_workspace * w)
{
  if (w->movstat_workspace_p)
    gsl_movstat_free(w->movstat_workspace_p);

  free(w);
}

/*
gsl_filter_median()
  Standard median filter

Inputs: endtype - end point handling
        x       - input vector
        y       - output vector
        w       - workspace
*/

int
gsl_filter_median(const gsl_filter_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_filter_median_workspace * w)
{
  int status = gsl_movstat_median((gsl_movstat_end_t) endtype, x, y, w->movstat_workspace_p);
  return status;
}
