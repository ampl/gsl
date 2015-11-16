/* rstat/rstat.c
 * 
 * Copyright (C) 2015 Patrick Alken
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
#include <gsl/gsl_rstat.h>

gsl_rstat_workspace *
gsl_rstat_alloc(void)
{
  gsl_rstat_workspace *w;

  w = calloc(1, sizeof(gsl_rstat_workspace));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->median_workspace_p = gsl_rstat_quantile_alloc(0.5);

  if (w->median_workspace_p == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for median workspace",
                      GSL_ENOMEM);
    }

  gsl_rstat_reset(w);

  return w;
} /* gsl_rstat_alloc() */

void
gsl_rstat_free(gsl_rstat_workspace *w)
{
  if (w->median_workspace_p)
    gsl_rstat_quantile_free(w->median_workspace_p);

  free(w);
} /* gsl_rstat_free() */

size_t
gsl_rstat_n(gsl_rstat_workspace *w)
{
  return w->n;
} /* gsl_rstat_n() */

/* add a data point to the running totals */
int
gsl_rstat_add(const double x, gsl_rstat_workspace *w)
{
  double delta = x - w->mean;
  double delta_n, delta_nsq, term1, n;

  /* update min and max */
  if (w->n == 0)
    {
      w->min = x;
      w->max = x;
    }
  else
    {
      if (x < w->min)
        w->min = x;
      if (x > w->max)
        w->max = x;
    }

  /* update mean and variance */
  n = (double) ++(w->n);
  delta_n = delta / n;
  delta_nsq = delta_n * delta_n;
  term1 = delta * delta_n * (n - 1.0);
  w->mean += delta_n;
  w->M4 += term1 * delta_nsq * (n * n - 3.0 * n + 3.0) +
           6.0 * delta_nsq * w->M2 - 4.0 * delta_n * w->M3;
  w->M3 += term1 * delta_n * (n - 2.0) - 3.0 * delta_n * w->M2;
  w->M2 += term1;

  /* update median */
  gsl_rstat_quantile_add(x, w->median_workspace_p);

  return GSL_SUCCESS;
} /* gsl_rstat_add() */

double
gsl_rstat_min(gsl_rstat_workspace *w)
{
  return w->min;
} /* gsl_rstat_min() */

double
gsl_rstat_max(gsl_rstat_workspace *w)
{
  return w->max;
} /* gsl_rstat_max() */

double
gsl_rstat_mean(gsl_rstat_workspace *w)
{
  return w->mean;
} /* gsl_rstat_mean() */

double
gsl_rstat_variance(gsl_rstat_workspace *w)
{
  if (w->n > 1)
    {
      double n = (double) w->n;
      return (w->M2 / (n - 1.0));
    }
  else
    return 0.0;
} /* gsl_rstat_variance() */

double
gsl_rstat_sd(gsl_rstat_workspace *w)
{
  double var = gsl_rstat_variance(w);

  return (sqrt(var));
} /* gsl_rstat_sd() */

/* standard deviation of the mean: sigma / sqrt(n) */
double
gsl_rstat_sd_mean(gsl_rstat_workspace *w)
{
  if (w->n > 0)
    {
      double sd = gsl_rstat_sd(w);
      return (sd / sqrt((double) w->n));
    }
  else
    return 0.0;
} /* gsl_rstat_sd_mean() */

double
gsl_rstat_median(gsl_rstat_workspace *w)
{
  return gsl_rstat_quantile_get(w->median_workspace_p);
}

double
gsl_rstat_skew(gsl_rstat_workspace *w)
{
  if (w->n > 0)
    {
      double n = (double) w->n;
      double fac = pow(n - 1.0, 1.5) / n;
      return ((fac * w->M3) / pow(w->M2, 1.5));
    }
  else
    return 0.0;
} /* gsl_rstat_skew() */

double
gsl_rstat_kurtosis(gsl_rstat_workspace *w)
{
  if (w->n > 0)
    {
      double n = (double) w->n;
      double fac = ((n - 1.0) / n) * (n - 1.0);
      return ((fac * w->M4) / (w->M2 * w->M2) - 3.0);
    }
  else
    return 0.0;
} /* gsl_rstat_kurtosis() */

int
gsl_rstat_reset(gsl_rstat_workspace *w)
{
  w->min = 0.0;
  w->max = 0.0;
  w->mean = 0.0;
  w->M2 = 0.0;
  w->M3 = 0.0;
  w->M4 = 0.0;
  w->n = 0;

  return GSL_SUCCESS;
} /* gsl_rstat_reset() */
