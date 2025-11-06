/* filter/rmedian.c
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

/*
 * This module contains routines related to the recursive median filter. The
 * algorithm is based on this paper,
 *
 * [1] S-J Ko, Y. H. Lee, and A. T. Fam, Efficient Implementation of One-Dimensional
 * Recursive Median Filters, IEEE Transactions on Circuits and Systems, Vol 37,
 * No 11, 1990.
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_movstat.h>

typedef struct
{
  const gsl_movstat_accum * minmax_acc; /* minimum/maximum accumulator */
  void *minmax_state;                   /* minimum/maximum accumulator workspace */
} rmedian_state_t;

static size_t rmedian_size(const size_t n);
static int rmedian_init(const size_t n, void * vstate);
static int rmedian_insert(const double x, void * vstate);
static int rmedian_delete(void * vstate);
static int rmedian_get(void * params, double * result, const void * vstate);

static const gsl_movstat_accum rmedian_accum_type;

gsl_filter_rmedian_workspace *
gsl_filter_rmedian_alloc(const size_t K)
{
  gsl_filter_rmedian_workspace *w;
  size_t state_size;

  w = calloc(1, sizeof(gsl_filter_rmedian_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->H = K / 2;
  w->K = 2*w->H + 1;
  w->minmaxacc = gsl_movstat_accum_minmax;

  w->window = malloc(w->K * sizeof(double));
  if (w->window == NULL)
    {
      gsl_filter_rmedian_free(w);
      GSL_ERROR_NULL ("failed to allocate space for window", GSL_ENOMEM);
    }

  state_size = rmedian_size(w->H + 1);

  w->state = malloc(state_size);
  if (w->state == NULL)
    {
      gsl_filter_rmedian_free(w);
      GSL_ERROR_NULL ("failed to allocate space for min/max state", GSL_ENOMEM);
    }

  w->movstat_workspace_p = gsl_movstat_alloc_with_size(state_size, 0, w->H);
  if (!w->movstat_workspace_p)
    {
      gsl_filter_rmedian_free(w);
      GSL_ERROR_NULL ("failed to allocate space for movstat workspace", GSL_ENOMEM);
    }

  return w;
}

void
gsl_filter_rmedian_free(gsl_filter_rmedian_workspace * w)
{
  if (w->state)
    free(w->state);

  if (w->window)
    free(w->window);

  if (w->movstat_workspace_p)
    gsl_movstat_free(w->movstat_workspace_p);

  free(w);
}

/*
gsl_filter_rmedian()
  Recursive median filter

Inputs: endtype - end point handling
        x       - input vector
        y       - output vector
        w       - workspace
*/

int
gsl_filter_rmedian(const gsl_filter_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_filter_rmedian_workspace * w)
{
  if (x->size != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      int status = GSL_SUCCESS;
      const size_t n = x->size;
      const int H = (int) w->H;
      double yprev;
      int wsize;

      /* find median of first window to initialize filter */
      wsize = gsl_movstat_fill((gsl_movstat_end_t) endtype, x, 0, H, H, w->window);
      yprev = gsl_stats_median(w->window, 1, wsize);
      gsl_vector_set(y, 0, yprev);

      if (x->size > 1)
        {
          gsl_vector_const_view xv = gsl_vector_const_subvector(x, 1, n - 1);
          gsl_vector_view yv = gsl_vector_subvector(y, 1, n - 1);

          /* apply recursive median filter to x[2:end] */
          status = gsl_movstat_apply_accum((gsl_movstat_end_t) endtype, &xv.vector, &rmedian_accum_type, (void *) &yprev, &yv.vector,
                                           NULL, w->movstat_workspace_p);
        }

      return status;
    }
}

static size_t
rmedian_size(const size_t n)
{
  size_t size = 0;
  const gsl_movstat_accum * acc = gsl_movstat_accum_minmax;

  size += sizeof(rmedian_state_t);
  size += (acc->size)(n);

  return size;
}

static int
rmedian_init(const size_t n, void * vstate)
{
  rmedian_state_t * state = (rmedian_state_t *) vstate;

  state->minmax_acc = gsl_movstat_accum_minmax;
  state->minmax_state = (void *) ((unsigned char *) vstate + sizeof(rmedian_state_t));

  (state->minmax_acc->init)(n, state->minmax_state);

  return GSL_SUCCESS;
}

static int
rmedian_insert(const double x, void * vstate)
{
  rmedian_state_t * state = (rmedian_state_t *) vstate;
  return (state->minmax_acc->insert)(x, state->minmax_state);
}

static int
rmedian_delete(void * vstate)
{
  rmedian_state_t * state = (rmedian_state_t *) vstate;
  return (state->minmax_acc->delete_oldest)(state->minmax_state);
}

static int
rmedian_get(void * params, double * result, const void * vstate)
{
  const rmedian_state_t * state = (const rmedian_state_t *) vstate;
  double *yprev = (double *) params; /* previous filter output */
  double y;                          /* new filter output */
  double xminmax[2];

  /* get minimum/maximum values of {x_i,...,x_{i+H}} */
  (state->minmax_acc->get)(NULL, xminmax, state->minmax_state);

  /* y = median [ yprev, xmin, xmax ] */
  if (*yprev <= xminmax[0])
    y = xminmax[0];
  else if (*yprev <= xminmax[1])
    y = *yprev;
  else
    y = xminmax[1];

  *result = y;
  *yprev = y;

  return GSL_SUCCESS;
}

static const gsl_movstat_accum rmedian_accum_type =
{
  rmedian_size,
  rmedian_init,
  rmedian_insert,
  rmedian_delete,
  rmedian_get
};
