/* movstat/mvacc.c
 *
 * Moving window mean/variance accumulator - based on a modification
 * to Welford's algorithm, discussed here:
 *
 * https://stackoverflow.com/a/6664212
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

typedef double ringbuf_type_t;

#include "ringbuf.c"

typedef struct
{
  size_t n;      /* window size */
  size_t k;      /* number of samples currently in window */
  double mean;   /* current window mean */
  double M2;     /* current window M2 */
  ringbuf *rbuf; /* ring buffer storing current window */
} mvacc_state_t;

static size_t
mvacc_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(mvacc_state_t);
  size += ringbuf_size(n);

  return size;
}

static int
mvacc_init(const size_t n, void * vstate)
{
  mvacc_state_t * state = (mvacc_state_t *) vstate;

  state->n = n;
  state->k = 0;
  state->mean = 0.0;
  state->M2 = 0.0;

  state->rbuf = (ringbuf *) ((unsigned char *) vstate + sizeof(mvacc_state_t));
  ringbuf_init(n, state->rbuf);

  return GSL_SUCCESS;
}

static int
mvacc_insert(const double x, void * vstate)
{
  mvacc_state_t * state = (mvacc_state_t *) vstate;

  if (ringbuf_is_full(state->rbuf))
    {
      /* remove oldest window element and add new one */
      double old = ringbuf_peek_back(state->rbuf);
      double prev_mean = state->mean;

      state->mean += (x - old) / (double) state->n;
      state->M2 += ((old - prev_mean) + (x - state->mean)) * (x - old);
    }
  else
    {
      double delta = x - state->mean;

      /*
       * Welford algorithm:
       *
       * mu_new = mu_old + (x - mu_old) / n
       * M2_new = M2_old + (x - mu_old) * (x - mu_new)
       */

      ++(state->k);
      state->mean += delta / (double) state->k;
      state->M2 += delta * (x - state->mean);
    }

  /* add new element to ring buffer */
  ringbuf_insert(x, state->rbuf);

  return GSL_SUCCESS;
}

static int
mvacc_delete(void * vstate)
{
  mvacc_state_t * state = (mvacc_state_t *) vstate;

  if (!ringbuf_is_empty(state->rbuf))
    {
      if (state->k > 1)
        {
          /*
           * mu_new = mu_old + (mu_old - x_old) / (n - 1)
           * M2_new = M2_old - (mu_old - x_old) * (mu_new - x_old)
           */

          double old = ringbuf_peek_back(state->rbuf);
          double prev_mean = state->mean;
          double delta = prev_mean - old;

          state->mean += delta / (state->k - 1.0);
          state->M2 -= delta * (state->mean - old);
        }
      else if (state->k == 1)
        {
          state->mean = 0.0;
          state->M2 = 0.0;
        }

      ringbuf_pop_back(state->rbuf);
      --(state->k);
    }

  return GSL_SUCCESS;
}

static int
mvacc_mean(void * params, double * result, const void * vstate)
{
  const mvacc_state_t * state = (const mvacc_state_t *) vstate;
  (void) params;
  *result = state->mean;
  return GSL_SUCCESS;
}

static int
mvacc_variance(void * params, double * result, const void * vstate)
{
  const mvacc_state_t * state = (const mvacc_state_t *) vstate;

  (void) params;

  if (state->k < 2)
    *result = 0.0;
  else
    *result = state->M2 / (state->k - 1.0);

  return GSL_SUCCESS;
}

static int
mvacc_sd(void * params, double * result, const void * vstate)
{
  double variance;
  int status = mvacc_variance(params, &variance, vstate);
  *result = sqrt(variance);
  return status;
}

static const gsl_movstat_accum mean_accum_type =
{
  mvacc_size,
  mvacc_init,
  mvacc_insert,
  mvacc_delete,
  mvacc_mean
};

const gsl_movstat_accum *gsl_movstat_accum_mean = &mean_accum_type;

static const gsl_movstat_accum variance_accum_type =
{
  mvacc_size,
  mvacc_init,
  mvacc_insert,
  mvacc_delete,
  mvacc_variance
};

const gsl_movstat_accum *gsl_movstat_accum_variance = &variance_accum_type;

static const gsl_movstat_accum sd_accum_type =
{
  mvacc_size,
  mvacc_init,
  mvacc_insert,
  mvacc_delete,
  mvacc_sd
};

const gsl_movstat_accum *gsl_movstat_accum_sd = &sd_accum_type;
