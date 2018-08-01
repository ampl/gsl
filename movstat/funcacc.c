/* movstat/funcacc.c
 *
 * Moving window accumulator for arbitrary user-defined function
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

typedef double funcacc_type_t;
typedef funcacc_type_t ringbuf_type_t;

#include "ringbuf.c"

typedef struct
{
  funcacc_type_t *window; /* linear array for current window */
  ringbuf *rbuf;          /* ring buffer storing current window */
} funcacc_state_t;

static size_t
funcacc_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(funcacc_state_t);
  size += n * sizeof(funcacc_type_t);
  size += ringbuf_size(n);

  return size;
}

static int
funcacc_init(const size_t n, void * vstate)
{
  funcacc_state_t * state = (funcacc_state_t *) vstate;

  state->window = (funcacc_type_t *) ((unsigned char *) vstate + sizeof(funcacc_state_t));
  state->rbuf = (ringbuf *) ((unsigned char *) state->window + n * sizeof(funcacc_type_t));

  ringbuf_init(n, state->rbuf);

  return GSL_SUCCESS;
}

static int
funcacc_insert(const funcacc_type_t x, void * vstate)
{
  funcacc_state_t * state = (funcacc_state_t *) vstate;

  /* add new element to ring buffer */
  ringbuf_insert(x, state->rbuf);

  return GSL_SUCCESS;
}

static int
funcacc_delete(void * vstate)
{
  funcacc_state_t * state = (funcacc_state_t *) vstate;

  if (!ringbuf_is_empty(state->rbuf))
    ringbuf_pop_back(state->rbuf);

  return GSL_SUCCESS;
}

static int
funcacc_get(void * params, funcacc_type_t * result, const void * vstate)
{
  const funcacc_state_t * state = (const funcacc_state_t *) vstate;
  gsl_movstat_function *f = (gsl_movstat_function *) params;
  size_t n = ringbuf_copy(state->window, state->rbuf);

  *result = GSL_MOVSTAT_FN_EVAL(f, n, state->window);

  return GSL_SUCCESS;
}

static const gsl_movstat_accum func_accum_type =
{
  funcacc_size,
  funcacc_init,
  funcacc_insert,
  funcacc_delete,
  funcacc_get
};

const gsl_movstat_accum *gsl_movstat_accum_userfunc = &func_accum_type;
