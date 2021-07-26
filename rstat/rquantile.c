/* rstat/rquantile.c
 * 
 * Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021 Patrick Alken
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
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>

/*
 * Running quantile calculation based on the paper
 *
 * [1] R. Jain and I. Chlamtac, "The P^2 algorithm for dynamic
 *     calculation of quantiles and histograms without storing
 *     observations", Communications of the ACM, October 1985
 */

static double calc_psq(const double qp1, const double q, const double qm1,
                       const double d, const double np1, const double n, const double nm1);

gsl_rstat_quantile_workspace *
gsl_rstat_quantile_alloc(const double p)
{
  gsl_rstat_quantile_workspace *w;

  w = calloc(1, sizeof(gsl_rstat_quantile_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->p = p;

  gsl_rstat_quantile_reset(w);

  return w;
} /* gsl_rstat_quantile_alloc() */

void
gsl_rstat_quantile_free(gsl_rstat_quantile_workspace *w)
{
  free(w);
} /* gsl_rstat_quantile_free() */

int
gsl_rstat_quantile_reset(gsl_rstat_quantile_workspace *w)
{
  const double p = w->p;
  size_t i;

  /* initialize positions n */
  for (i = 0; i < 5; ++i)
    w->npos[i] = i + 1;

  /* initialize n' */
  w->np[0] = 1.0;
  w->np[1] = 1.0 + 2.0 * p;
  w->np[2] = 1.0 + 4.0 * p;
  w->np[3] = 3.0 + 2.0 * p;
  w->np[4] = 5.0;

  /* initialize dn' */
  w->dnp[0] = 0.0;
  w->dnp[1] = 0.5 * p;
  w->dnp[2] = p;
  w->dnp[3] = 0.5 * (1.0 + p);
  w->dnp[4] = 1.0;

  w->n = 0;

  return GSL_SUCCESS;
}

int
gsl_rstat_quantile_add(const double x, gsl_rstat_quantile_workspace *w)
{
  if (w->n < 5)
    {
      w->q[w->n] = x;
    }
  else
    {
      int i;
      int k = -1;

      if (w->n == 5)
        {
          /* initialization: sort the first five heights */
          gsl_sort(w->q, 1, w->n);
        }

      /* step B1: find k such that q_k <= x < q_{k+1} */
      if (x < w->q[0])
        {
          w->q[0] = x;
          k = 0;
        }
      else if (x >= w->q[4])
        {
          w->q[4] = x;
          k = 3;
        }
      else
        {
          for (i = 0; i <= 3; ++i)
            {
              if (w->q[i] <= x && x < w->q[i + 1])
                {
                  k = i;
                  break;
                }
            }
        }

      if (k < 0)
        {
          /* we could get here if x is nan */
          GSL_ERROR ("invalid input argument x", GSL_EINVAL);
        }

      /* step B2(a): update n_i */
      for (i = k + 1; i <= 4; ++i)
        ++(w->npos[i]);

      /* step B2(b): update n_i' */
      for (i = 0; i < 5; ++i)
        w->np[i] += w->dnp[i];

      /* step B3: update heights */
      for (i = 1; i <= 3; ++i)
        {
          double ni = (double) w->npos[i];
          double d = w->np[i] - ni;

          if ((d >= 1.0 && (w->npos[i + 1] - w->npos[i] > 1)) ||
              (d <= -1.0 && (w->npos[i - 1] - w->npos[i] < -1)))
            {
              int dsign = (d > 0.0) ? 1 : -1;
              double qp1 = w->q[i + 1];
              double qi = w->q[i];
              double qm1 = w->q[i - 1];
              double np1 = (double) w->npos[i + 1];
              double nm1 = (double) w->npos[i - 1];
              double qp = calc_psq(qp1, qi, qm1, (double) dsign,
                                   np1, ni, nm1);

              if (qm1 < qp && qp < qp1)
                w->q[i] = qp;
              else
                {
                  /* use linear formula */
                  w->q[i] += dsign * (w->q[i + dsign] - qi) / ((double) w->npos[i + dsign] - ni);
                }

              w->npos[i] += dsign;
            }
        }
    }

  ++(w->n);

  return GSL_SUCCESS;
} /* gsl_rstat_quantile_add() */

double
gsl_rstat_quantile_get(gsl_rstat_quantile_workspace *w)
{
  if (w->n > 5)
    {
      return w->q[2];
    }
  else
    {
      /* not yet initialized */
      gsl_sort(w->q, 1, w->n);
      return gsl_stats_quantile_from_sorted_data(w->q, 1, w->n, w->p);
    }
} /* gsl_rstat_quantile_get() */

static double
calc_psq(const double qp1, const double q, const double qm1,
         const double d, const double np1, const double n, const double nm1)
{
  double outer = d / (np1 - nm1);
  double inner_left = (n - nm1 + d) * (qp1 - q) / (np1 - n);
  double inner_right = (np1 - n - d) * (q - qm1) / (n - nm1);

  return q + outer * (inner_left + inner_right);
} /* calc_psq() */
