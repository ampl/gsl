/* filter/impulse.c
 *
 * Impulse detecting filters
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

static int filter_impulse(const double scale, const double epsilon, const double t, const gsl_vector * x, const gsl_vector * xmedian,
                          gsl_vector * y, gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier);
 
/*
gsl_filter_impulse_alloc()
  Allocate a workspace for impulse detection filtering.

Inputs: K - number of samples in window; if even, it is rounded up to
            the next odd, to have a symmetric window

Return: pointer to workspace
*/

gsl_filter_impulse_workspace *
gsl_filter_impulse_alloc(const size_t K)
{
  gsl_filter_impulse_workspace *w;

  w = calloc(1, sizeof(gsl_filter_impulse_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->movstat_workspace_p = gsl_movstat_alloc(K);
  if (w->movstat_workspace_p == 0)
    {
      gsl_filter_impulse_free(w);
      return NULL;
    }

  return w;
}

void
gsl_filter_impulse_free(gsl_filter_impulse_workspace * w)
{
  if (w->movstat_workspace_p)
    gsl_movstat_free(w->movstat_workspace_p);

  free(w);
}

/*
gsl_filter_impulse()
  Apply an impulse detection filter to an input vector. The filter output is

y_i = { x_i, |x_i - m_i| <= t * S_i
      { m_i, |x_i - m_i| > t * S_i

where m_i is the median of the window W_i^H and S_i is the scale estimate (MAD, IQR, S_n, Q_n)

Inputs: endtype    - how to handle signal end points
        scale_type - which statistic to use for scale estimate (MAD, IQR, etc)
        t          - number of standard deviations required to identity outliers (>= 0)
        x          - input vector, size n
        y          - (output) filtered vector, size n
        xmedian    - (output) vector of median values of x, size n
                     xmedian_i = median of window centered on x_i
        xsigma     - (output) vector of estimated local standard deviations of x, size n
                     xsigma_i = sigma for i-th window: scale*MAD
        noutlier   - (output) number of outliers detected
        ioutlier   - (output) boolean array indicating outliers identified, size n; may be NULL
                     ioutlier_i = 1 if outlier detected, 0 if not
        w          - workspace

Notes:
*/

int
gsl_filter_impulse(const gsl_filter_end_t endtype, const gsl_filter_scale_t scale_type, const double t,
                   const gsl_vector * x, gsl_vector * y, gsl_vector * xmedian, gsl_vector * xsigma, size_t * noutlier,
                   gsl_vector_int * ioutlier, gsl_filter_impulse_workspace * w)
{
  const size_t n = x->size;

  if (n != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else if (xmedian->size != n)
    {
      GSL_ERROR("xmedian vector must match input size", GSL_EBADLEN);
    }
  else if (xsigma->size != n)
    {
      GSL_ERROR("xsigma vector must match input size", GSL_EBADLEN);
    }
  else if ((ioutlier != NULL) && (ioutlier->size != n))
    {
      GSL_ERROR("ioutlier vector must match input size", GSL_EBADLEN);
    }
  else if (t < 0.0)
    {
      GSL_ERROR("t must be non-negative", GSL_EDOM);
    }
  else
    {
      int status;
      double scale = 1.0;

      switch (scale_type)
        {
          case GSL_FILTER_SCALE_MAD:
            {
              /* compute window medians and MADs */
              gsl_movstat_mad((gsl_movstat_end_t) endtype, x, xmedian, xsigma, w->movstat_workspace_p);

              break;
            }

          case GSL_FILTER_SCALE_IQR:
            {
              /* multiplication factor for IQR to estimate stddev for Gaussian signal */
              scale = 0.741301109252801;

              /* calculate the window medians */
              gsl_movstat_median((gsl_movstat_end_t) endtype, x, xmedian, w->movstat_workspace_p);
      
              /* calculate window IQRs */
              gsl_movstat_qqr((gsl_movstat_end_t) endtype, x, 0.25, xsigma, w->movstat_workspace_p);

              break;
            }

          case GSL_FILTER_SCALE_SN:
            {
              /* calculate the window medians */
              gsl_movstat_median((gsl_movstat_end_t) endtype, x, xmedian, w->movstat_workspace_p);
      
              /* calculate window S_n values */
              gsl_movstat_Sn((gsl_movstat_end_t) endtype, x, xsigma, w->movstat_workspace_p);

              break;
            }

          case GSL_FILTER_SCALE_QN:
            {
              /* calculate the window medians */
              gsl_movstat_median((gsl_movstat_end_t) endtype, x, xmedian, w->movstat_workspace_p);
      
              /* calculate window Q_n values */
              gsl_movstat_Qn((gsl_movstat_end_t) endtype, x, xsigma, w->movstat_workspace_p);

              break;
            }

          default:
            GSL_ERROR("unknown scale type", GSL_EDOM);
            break;
        }

      /* apply impulse detecting filter using previously computed scale estimate */
      status = filter_impulse(scale, 0.0, t, x, xmedian, y, xsigma, noutlier, ioutlier);

      return status;
    }
}

/*
filter_impulse()
  Apply an impulse detection filter to an input vector. The filter output is

y_i = { x_i, |x_i - m_i| <= t * S_i OR S_i < epsilon
      { m_i, |x_i - m_i| > t * S_i

where m_i is the median of the window W_i^H and S_i is the scale estimate (MAD, IQR, etc)

Inputs: scale    - scale factor to multiply xsigma to get unbiased estimate of stddev for Gaussian data
        epsilon  - minimum allowed scale estimate for identifying outliers
        t        - number of standard deviations required to identity outliers (>= 0)
        x        - input vector, size n
        xmedian  - vector of median values of x, size n
                   xmedian_i = median of window centered on x_i
        y        - (output) filtered vector, size n
        xsigma   - (output) vector of estimated local standard deviations of x, size n
                   xsigma_i = S_n for i-th window
        noutlier - (output) number of outliers detected
        ioutlier - (output) boolean array indicating outliers identified, size n; may be NULL
                   ioutlier_i = 1 if outlier detected, 0 if not

Notes:
1) If S_i = 0 or is very small for a particular sample, then the filter may erroneously flag the
sample as an outlier, since it will act as a standard median filter. To avoid this scenario, the
parameter epsilon specifies the minimum value of S_i which can be used in the filter test. Any
samples for which S_i < epsilon are passed through unchanged.
*/

static int
filter_impulse(const double scale, const double epsilon, const double t, const gsl_vector * x, const gsl_vector * xmedian,
               gsl_vector * y, gsl_vector * xsigma, size_t * noutlier, gsl_vector_int * ioutlier)
{
  const size_t n = x->size;

  if (n != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else if (xmedian->size != n)
    {
      GSL_ERROR("xmedian vector must match input size", GSL_EBADLEN);
    }
  else if (xsigma->size != n)
    {
      GSL_ERROR("xsigma vector must match input size", GSL_EBADLEN);
    }
  else if ((ioutlier != NULL) && (ioutlier->size != n))
    {
      GSL_ERROR("ioutlier vector must match input size", GSL_EBADLEN);
    }
  else if (t < 0.0)
    {
      GSL_ERROR("t must be non-negative", GSL_EDOM);
    }
  else
    {
      size_t i;

      *noutlier = 0;

      /* build output vector */
      for (i = 0; i < n; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double xmedi = gsl_vector_get(xmedian, i);
          double absdevi = fabs(xi - xmedi); /* absolute deviation for this sample */
          double *xsigmai = gsl_vector_ptr(xsigma, i);

          /* multiply by scale factor to get estimate of standard deviation */
          *xsigmai *= scale;

          /*
           * If the absolute deviation for this sample is more than t stddevs
           * for this window (and S_i is sufficiently large to avoid scale implosion),
           * set the output value to the window median, otherwise use the original sample
           */
          if ((*xsigmai >= epsilon) && (absdevi > t * (*xsigmai)))
            {
              gsl_vector_set(y, i, xmedi);
              ++(*noutlier);
              if (ioutlier)
                gsl_vector_int_set(ioutlier, i, 1);
            }
          else
            {
              gsl_vector_set(y, i, xi);
              if (ioutlier)
                gsl_vector_int_set(ioutlier, i, 0);
            }
        }

      return GSL_SUCCESS;
    }
}
