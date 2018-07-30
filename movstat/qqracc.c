/* movstat/qqracc.c
 *
 * Moving window QQR accumulator
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

typedef double qqracc_type_t;
typedef qqracc_type_t ringbuf_type_t;

#include "ringbuf.c"

typedef struct
{
  qqracc_type_t *window; /* linear array for current window */
  ringbuf *rbuf;         /* ring buffer storing current window */
} qqracc_state_t;

static size_t
qqracc_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(qqracc_state_t);
  size += n * sizeof(qqracc_type_t);
  size += ringbuf_size(n);

  return size;
}

static int
qqracc_init(const size_t n, void * vstate)
{
  qqracc_state_t * state = (qqracc_state_t *) vstate;

  state->window = (qqracc_type_t *) ((unsigned char *) vstate + sizeof(qqracc_state_t));
  state->rbuf = (ringbuf *) ((unsigned char *) state->window + n * sizeof(qqracc_type_t));

  ringbuf_init(n, state->rbuf);

  return GSL_SUCCESS;
}

static int
qqracc_insert(const qqracc_type_t x, void * vstate)
{
  qqracc_state_t * state = (qqracc_state_t *) vstate;

  /* add new element to ring buffer */
  ringbuf_insert(x, state->rbuf);

  return GSL_SUCCESS;
}

static int
qqracc_delete(void * vstate)
{
  qqracc_state_t * state = (qqracc_state_t *) vstate;

  if (!ringbuf_is_empty(state->rbuf))
    ringbuf_pop_back(state->rbuf);

  return GSL_SUCCESS;
}

/* FIXME XXX: this is inefficient - could be improved by maintaining a sorted ring buffer */
static int
qqracc_get(void * params, qqracc_type_t * result, const void * vstate)
{
  const qqracc_state_t * state = (const qqracc_state_t *) vstate;
  double q = *(double *) params;
  size_t n = ringbuf_copy(state->window, state->rbuf);
  double quant1, quant2;

  gsl_sort(state->window, 1, n);

  /* compute q-quantile and (1-q)-quantile */
  quant1 = gsl_stats_quantile_from_sorted_data(state->window, 1, n, q);
  quant2 = gsl_stats_quantile_from_sorted_data(state->window, 1, n, 1.0 - q);

  /* compute q-quantile range */
  *result = quant2 - quant1;

  return GSL_SUCCESS;
}

static const gsl_movstat_accum qqr_accum_type =
{
  qqracc_size,
  qqracc_init,
  qqracc_insert,
  qqracc_delete,
  qqracc_get
};

const gsl_movstat_accum *gsl_movstat_accum_qqr = &qqr_accum_type;
