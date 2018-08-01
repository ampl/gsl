/* movstat/snacc.c
 *
 * Moving window S_n accumulator
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

typedef double snacc_type_t;
typedef snacc_type_t ringbuf_type_t;

#include "ringbuf.c"

typedef struct
{
  snacc_type_t *window; /* linear array for current window */
  snacc_type_t *work;   /* workspace */
  ringbuf *rbuf;        /* ring buffer storing current window */
} snacc_state_t;

static size_t
snacc_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(snacc_state_t);
  size += 2 * n * sizeof(snacc_type_t);
  size += ringbuf_size(n);

  return size;
}

static int
snacc_init(const size_t n, void * vstate)
{
  snacc_state_t * state = (snacc_state_t *) vstate;

  state->window = (snacc_type_t *) ((unsigned char *) vstate + sizeof(snacc_state_t));
  state->work = (snacc_type_t *) ((unsigned char *) state->window + n * sizeof(snacc_type_t));
  state->rbuf = (ringbuf *) ((unsigned char *) state->work + n * sizeof(snacc_type_t));

  ringbuf_init(n, state->rbuf);

  return GSL_SUCCESS;
}

static int
snacc_insert(const snacc_type_t x, void * vstate)
{
  snacc_state_t * state = (snacc_state_t *) vstate;

  /* add new element to ring buffer */
  ringbuf_insert(x, state->rbuf);

  return GSL_SUCCESS;
}

static int
snacc_delete(void * vstate)
{
  snacc_state_t * state = (snacc_state_t *) vstate;

  if (!ringbuf_is_empty(state->rbuf))
    ringbuf_pop_back(state->rbuf);

  return GSL_SUCCESS;
}

/* FIXME XXX: this is inefficient - could be improved by maintaining a sorted ring buffer */
static int
snacc_get(void * params, snacc_type_t * result, const void * vstate)
{
  const snacc_state_t * state = (const snacc_state_t *) vstate;
  size_t n = ringbuf_copy(state->window, state->rbuf);

  (void) params;

  gsl_sort(state->window, 1, n);
  *result = gsl_stats_Sn_from_sorted_data(state->window, 1, n, state->work);

  return GSL_SUCCESS;
}

static const gsl_movstat_accum sn_accum_type =
{
  snacc_size,
  snacc_init,
  snacc_insert,
  snacc_delete,
  snacc_get
};

const gsl_movstat_accum *gsl_movstat_accum_Sn = &sn_accum_type;
