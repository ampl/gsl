/* movstat/mmacc.c
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
 * This module contains routines for tracking minimum/maximum values of a
 * moving fixed-sized window. It is based on the algorithm of:
 *
 * [1] Daniel Lemire, Streaming Maximum-Minimum Filter Using No More than Three Comparisons per Element,
 *     Nordic Journal of Computing, Volume 13, Number 4, pages 328-339, 2006
 *
 * Also available as a preprint here: https://arxiv.org/abs/cs/0610046
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>

typedef double mmacc_type_t;
typedef mmacc_type_t ringbuf_type_t;

#include "deque.c"
#include "ringbuf.c"

typedef struct
{
  size_t n;             /* window size */
  size_t k;             /* number of samples in current window */
  mmacc_type_t xprev;   /* previous sample added to window */
  ringbuf *rbuf;        /* ring buffer storing current window, size n */
  deque *minque;        /* double-ended queue of min values (L) */
  deque *maxque;        /* double-ended queue of max values (U) */
} mmacc_state_t;

static size_t mmacc_size(const size_t n);
static int mmacc_init(const size_t n, void * vstate);
static int mmacc_insert(const mmacc_type_t x, void * vstate);
static int mmacc_delete(void * vstate);
static int mmacc_min(void * params, mmacc_type_t * result, const void * vstate);
static int mmacc_max(void * params, mmacc_type_t * result, const void * vstate);
static int mmacc_minmax(void * params, mmacc_type_t * result, const void * vstate);

static size_t
mmacc_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(mmacc_state_t);
  size += ringbuf_size(n);          /* rbuf */
  size += 2 * deque_size(n + 1);    /* minque/maxque */

  return size;
}

static int
mmacc_init(const size_t n, void * vstate)
{
  mmacc_state_t * state = (mmacc_state_t *) vstate;

  state->n = n;
  state->k = 0;
  state->xprev = 0.0;

  state->rbuf = (ringbuf *) ((unsigned char *) vstate + sizeof(mmacc_state_t));
  state->minque = (deque *) ((unsigned char *) state->rbuf + ringbuf_size(n));
  state->maxque = (deque *) ((unsigned char *) state->minque + deque_size(n + 1));

  ringbuf_init(n, state->rbuf);
  deque_init(n + 1, state->minque);
  deque_init(n + 1, state->maxque);

  return GSL_SUCCESS;
}

static int
mmacc_insert(const mmacc_type_t x, void * vstate)
{
  mmacc_state_t * state = (mmacc_state_t *) vstate;
  int head, tail;

  if (state->k == 0)
    {
      /* first sample */
      ringbuf_insert(x, state->rbuf);
      head = state->rbuf->head;
      deque_push_back(head, state->maxque);
      deque_push_back(head, state->minque);
    }
  else
    {
      if (x > state->xprev)
        {
          deque_pop_back(state->maxque);

          while (!deque_is_empty(state->maxque))
            {
              if (x <= state->rbuf->array[deque_peek_back(state->maxque)])
                break;

              deque_pop_back(state->maxque);
            }
        }
      else
        {
          deque_pop_back(state->minque);

          while (!deque_is_empty(state->minque))
            {
              if (x >= state->rbuf->array[deque_peek_back(state->minque)])
                break;

              deque_pop_back(state->minque);
            }
        }

      /* store new sample into ring buffer */
      tail = state->rbuf->tail;
      ringbuf_insert(x, state->rbuf);
      head = state->rbuf->head;

      deque_push_back(head, state->maxque);
      deque_push_back(head, state->minque);

      if (state->k == state->n)
        {
          /*
           * window is full - check if oldest window element is a global minimum/maximum
           * of current window - if so pop it from U/L queues;
           * the check head != tail ensures there is more than 1 element in the
           * queue, do not pop if queue has only 1 element, since this element would
           * be the newest sample
           */
          if (state->maxque->head != state->maxque->tail && tail == deque_peek_front(state->maxque))
            deque_pop_front(state->maxque);
          else if (state->minque->head != state->minque->tail && tail == deque_peek_front(state->minque))
            deque_pop_front(state->minque);
        }
    }

  if (state->k < state->n)
    ++(state->k);

  state->xprev = x;

  return GSL_SUCCESS;
}

static int
mmacc_delete(void * vstate)
{
  mmacc_state_t * state = (mmacc_state_t *) vstate;

  if (state->k > 0)
    {
      /*
       * check if oldest window element is a global minimum/maximum; if so
       * pop it from U/L queues
       */
      if (state->rbuf->tail == deque_peek_front(state->maxque))
        deque_pop_front(state->maxque);
      else if (state->rbuf->tail == deque_peek_front(state->minque))
        deque_pop_front(state->minque);

      /* remove oldest element from ring buffer */
      ringbuf_pop_back(state->rbuf);

      --(state->k);
    }

  return GSL_SUCCESS;
}

static int
mmacc_min(void * params, mmacc_type_t * result, const void * vstate)
{
  const mmacc_state_t * state = (const mmacc_state_t *) vstate;

  (void) params;

  if (state->k == 0)
    {
      GSL_ERROR ("no samples yet added to workspace", GSL_EINVAL);
    }
  else
    {
      *result = state->rbuf->array[deque_peek_front(state->minque)];
      return GSL_SUCCESS;
    }
}

static int
mmacc_max(void * params, mmacc_type_t * result, const void * vstate)
{
  const mmacc_state_t * state = (const mmacc_state_t *) vstate;

  (void) params;

  if (state->k == 0)
    {
      GSL_ERROR ("no samples yet added to workspace", GSL_EINVAL);
    }
  else
    {
      *result = state->rbuf->array[deque_peek_front(state->maxque)];
      return GSL_SUCCESS;
    }
}

static int
mmacc_minmax(void * params, mmacc_type_t * result, const void * vstate)
{
  const mmacc_state_t * state = (const mmacc_state_t *) vstate;

  (void) params;

  if (state->k == 0)
    {
      GSL_ERROR ("no samples yet added to workspace", GSL_EINVAL);
    }
  else
    {
      result[0] = state->rbuf->array[deque_peek_front(state->minque)];
      result[1] = state->rbuf->array[deque_peek_front(state->maxque)];
      return GSL_SUCCESS;
    }
}

static const gsl_movstat_accum min_accum_type =
{
  mmacc_size,
  mmacc_init,
  mmacc_insert,
  mmacc_delete,
  mmacc_min
};

const gsl_movstat_accum *gsl_movstat_accum_min = &min_accum_type;

static const gsl_movstat_accum max_accum_type =
{
  mmacc_size,
  mmacc_init,
  mmacc_insert,
  mmacc_delete,
  mmacc_max
};

const gsl_movstat_accum *gsl_movstat_accum_max = &max_accum_type;

static const gsl_movstat_accum minmax_accum_type =
{
  mmacc_size,
  mmacc_init,
  mmacc_insert,
  mmacc_delete,
  mmacc_minmax
};

const gsl_movstat_accum *gsl_movstat_accum_minmax = &minmax_accum_type;
