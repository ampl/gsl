/* movstat/madacc.c
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
 * This module contains routines for computing the median absolute deviation (MAD)
 * of a moving fixed-sized window.
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_statistics.h>

typedef double madacc_type_t;
typedef madacc_type_t ringbuf_type_t;

#include "ringbuf.c"

typedef struct
{
  size_t n;                        /* window size */
  const gsl_movstat_accum *medacc; /* median accumulator */
  void *median_state;              /* median accumulator workspace */
  ringbuf *rbuf;                   /* ring buffer storing current window, size n */
  madacc_type_t *work;             /* workspace, size n */
} madacc_state_t;

static size_t madacc_size(const size_t n);
static int madacc_init(const size_t n, void * vstate);
static int madacc_insert(const madacc_type_t x, void * vstate);
static int madacc_delete(void * vstate);
static int madacc_medmad(void * params, madacc_type_t * result, const void * vstate);

static size_t
madacc_size(const size_t n)
{
  size_t size = 0;
  const gsl_movstat_accum * acc = gsl_movstat_accum_median;

  size += sizeof(madacc_state_t);
  size += (acc->size)(n);            /* median accumulator */
  size += ringbuf_size(n);           /* rbuf */
  size += n * sizeof(madacc_type_t); /* work */

  return size;
}

static int
madacc_init(const size_t n, void * vstate)
{
  madacc_state_t * state = (madacc_state_t *) vstate;

  state->n = n;
  state->medacc = gsl_movstat_accum_median;

  state->median_state = (void *) ((unsigned char *) vstate + sizeof(madacc_state_t));
  state->rbuf = (ringbuf *) ((unsigned char *) state->median_state + (state->medacc->size)(n));
  state->work = (madacc_type_t *) ((unsigned char *) state->rbuf + ringbuf_size(n));

  (state->medacc->init)(n, state->median_state);
  ringbuf_init(n, state->rbuf);

  return GSL_SUCCESS;
}

static int
madacc_insert(const madacc_type_t x, void * vstate)
{
  madacc_state_t * state = (madacc_state_t *) vstate;

  /* insert element into median accumulator */
  (state->medacc->insert)(x, state->median_state);

  /* insert element into ring buffer */
  ringbuf_insert(x, state->rbuf);

  return GSL_SUCCESS;
}

static int
madacc_delete(void * vstate)
{
  madacc_state_t * state = (madacc_state_t *) vstate;

  if (!ringbuf_is_empty(state->rbuf))
    ringbuf_pop_back(state->rbuf);

  return GSL_SUCCESS;
}

static int
madacc_medmad(void * params, madacc_type_t * result, const void * vstate)
{
  const madacc_state_t * state = (const madacc_state_t *) vstate;

  if (ringbuf_is_empty(state->rbuf))
    {
      GSL_ERROR("no samples yet added to workspace", GSL_EINVAL);
    }
  else
    {
      int status;
      const double scale = *(double *) params;
      const int n = ringbuf_n(state->rbuf);
      int i;
      double median, mad;

      /* compute median of current window */
      status = (state->medacc->get)(NULL, &median, state->median_state);
      if (status)
        return status;

      /* compute current window absolute deviations from median */
      for (i = 0; i < n; ++i)
        {
          double xi = state->rbuf->array[(state->rbuf->head + i) % state->rbuf->size];
          state->work[i] = fabs(xi - median);
        }

      /* compute MAD of current window */
      mad = gsl_stats_median(state->work, 1, n);

      result[0] = median;
      result[1] = scale * mad;

      return GSL_SUCCESS;
    }
}

static const gsl_movstat_accum mad_accum_type =
{
  madacc_size,
  madacc_init,
  madacc_insert,
  NULL, /*XXX*/
  madacc_medmad
};

const gsl_movstat_accum *gsl_movstat_accum_mad = &mad_accum_type;
