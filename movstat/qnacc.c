/* movstat/qnacc.c
 *
 * Moving window Q_n accumulator
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

typedef double qnacc_type_t;
typedef qnacc_type_t ringbuf_type_t;

#include "ringbuf.c"

typedef struct
{
  qnacc_type_t *window; /* linear array for current window */
  qnacc_type_t *work;   /* workspace, length 3*n */
  int *work_int;        /* integer workspace, length 5*n */
  ringbuf *rbuf;        /* ring buffer storing current window */
} qnacc_state_t;

static size_t
qnacc_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(qnacc_state_t);
  size += n * sizeof(qnacc_type_t);     /* window */
  size += 3 * n * sizeof(qnacc_type_t); /* work */
  size += 5 * n * sizeof(int);          /* work_int */
  size += ringbuf_size(n);

  return size;
}

static int
qnacc_init(const size_t n, void * vstate)
{
  qnacc_state_t * state = (qnacc_state_t *) vstate;

  state->window = (qnacc_type_t *) ((unsigned char *) vstate + sizeof(qnacc_state_t));
  state->work = (qnacc_type_t *) ((unsigned char *) state->window + n * sizeof(qnacc_type_t));
  state->work_int = (int *) ((unsigned char *) state->work + 3 * n * sizeof(qnacc_type_t));
  state->rbuf = (ringbuf *) ((unsigned char *) state->work_int + 5 * n * sizeof(int));

  ringbuf_init(n, state->rbuf);

  return GSL_SUCCESS;
}

static int
qnacc_insert(const qnacc_type_t x, void * vstate)
{
  qnacc_state_t * state = (qnacc_state_t *) vstate;

  /* add new element to ring buffer */
  ringbuf_insert(x, state->rbuf);

  return GSL_SUCCESS;
}

static int
qnacc_delete(void * vstate)
{
  qnacc_state_t * state = (qnacc_state_t *) vstate;

  if (!ringbuf_is_empty(state->rbuf))
    ringbuf_pop_back(state->rbuf);

  return GSL_SUCCESS;
}

/* FIXME XXX: this is inefficient - could be improved by maintaining a sorted ring buffer */
static int
qnacc_get(void * params, qnacc_type_t * result, const void * vstate)
{
  const qnacc_state_t * state = (const qnacc_state_t *) vstate;
  size_t n = ringbuf_copy(state->window, state->rbuf);

  (void) params;

  gsl_sort(state->window, 1, n);
  *result = gsl_stats_Qn_from_sorted_data(state->window, 1, n, state->work, state->work_int);

  return GSL_SUCCESS;
}

static const gsl_movstat_accum qn_accum_type =
{
  qnacc_size,
  qnacc_init,
  qnacc_insert,
  qnacc_delete,
  qnacc_get
};

const gsl_movstat_accum *gsl_movstat_accum_Qn = &qn_accum_type;
