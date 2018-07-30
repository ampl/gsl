/* movstat/medacc.c
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
 *
 * Original copyright notice:
 * Copyright (c) 2011 ashelly.myopenid.com under <http://www.opensource.org/licenses/mit-license>
 */
 
#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_movstat.h>
 
#define ItemLess(a,b)  ((a)<(b))
#define ItemMean(a,b)  (((a)+(b))/2)
 
#define minCt(m) (((m)->ct-1)/2) /* count of items in minheap */
#define maxCt(m) (((m)->ct)/2)   /* count of items in maxheap */

typedef double medacc_type_t;

typedef struct
{
  int n;                /* window size */
  int idx;              /* position in circular queue */
  int ct;               /* count of items in queue */
  medacc_type_t *data;  /* circular queue of values, size k */
  int *pos;             /* index into `heap` for each value, size 2*k */
  int *heap;            /* max/median/min heap holding indices into `data` */
} medacc_state_t;

static size_t medacc_size(const size_t n);
static int medacc_init(const size_t n, void * vstate);
static int medacc_insert(const medacc_type_t x, void * vstate);
static int medacc_delete(void * vstate);
static int medacc_get(void * params, medacc_type_t * result, const void * vstate);

static int mmless(const medacc_state_t * state, const int i, const int j);
static int mmexchange(medacc_state_t * state, const int i, const int j);
static int mmCmpExch(medacc_state_t * state, const int i, const int j);
static void minSortDown(medacc_state_t * state, int i);
static void maxSortDown(medacc_state_t * state, int i);
static int minSortUp(medacc_state_t * state, int i);
static int maxSortUp(medacc_state_t * state, int i);

static size_t
medacc_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(medacc_state_t);
  size += n * sizeof(medacc_type_t);
  size += 2 * n * sizeof(int);

  return size;
}

static int
medacc_init(const size_t n, void * vstate)
{
  medacc_state_t * state = (medacc_state_t *) vstate;
  int k = (int) n;

  state->n = n;
  state->ct = 0;
  state->idx = 0;

  state->data = (medacc_type_t *) ((unsigned char *) vstate + sizeof(medacc_state_t));
  state->pos = (int *) ((unsigned char *) state->data + n * sizeof(medacc_type_t));
  state->heap = state->pos + n + (n/2); /* points to middle of storage */

  /* set up initial heap fill pattern: median,max,min,max,... */
  while (k--)
    {
      state->pos[k] = ((k + 1)/2) * ((k & 1) ? -1 : 1);
      state->heap[state->pos[k]] = k;
    }

  return GSL_SUCCESS;
}

static int
medacc_insert(const medacc_type_t x, void * vstate)
{
  medacc_state_t * state = (medacc_state_t *) vstate;
  int isNew = (state->ct < (int) state->n);
  int p = state->pos[state->idx];
  medacc_type_t old = state->data[state->idx];

  state->data[state->idx] = x;
  state->idx = (state->idx + 1) % state->n;
  state->ct += isNew;

  if (p > 0)       /* new item is in minHeap */
    {
      if (!isNew && ItemLess(old, x))
        minSortDown(state, p * 2);
      else if (minSortUp(state, p))
        maxSortDown(state, -1);
    }
  else if (p < 0)  /* new item is in maxHeap */
    {
      if (!isNew && ItemLess(x, old))
        maxSortDown(state, p * 2);
      else if (maxSortUp(state, p))
        minSortDown(state, 1);
    }
  else             /* new item is at median */
    {
      if (maxCt(state))
        maxSortDown(state, -1);

      if (minCt(state))
        minSortDown(state, 1);
    }

  return GSL_SUCCESS;
}

static int
medacc_delete(void * vstate)
{
  medacc_state_t * state = (medacc_state_t *) vstate;

  if (state->ct > 0)
    {
      int p = state->pos[(state->idx - state->ct + state->n) % state->n];

      if (p > 0)        /* oldest item is in minHeap */
        {
          mmexchange(state, p, minCt(state));
          --(state->ct);
          minSortDown(state, 2 * p);
        }
      else if (p < 0)   /* oldest item is in maxHeap */
        {
          mmexchange(state, p, maxCt(state));
          --(state->ct);
          maxSortDown(state, 2 * p);
        }
      else if (p == 0)  /* oldest item is at median */
        {
        }
    }

  return GSL_SUCCESS;
}

/* returns median (or average of 2 when item count is even) */
static int
medacc_get(void * params, medacc_type_t * result, const void * vstate)
{
  const medacc_state_t * state = (const medacc_state_t *) vstate;
  medacc_type_t median = state->data[state->heap[0]];

  (void) params;

  if ((state->ct & 1) == 0)
    median = ItemMean(median, state->data[state->heap[-1]]);

  *result = median;

  return GSL_SUCCESS;
}

/* returns 1 if heap[i] < heap[j] */
static int
mmless(const medacc_state_t * state, const int i, const int j)
{
  return ItemLess(state->data[state->heap[i]], state->data[state->heap[j]]);
}
 
/* swaps items i and j in heap, maintains indexes */
static int
mmexchange(medacc_state_t * state, const int i, const int j)
{
  int t = state->heap[i];
  state->heap[i] = state->heap[j];
  state->heap[j] = t;
  state->pos[state->heap[i]] = i;
  state->pos[state->heap[j]] = j;
  return 1;
}
 
/* swaps items i and j if i < j; returns true if swapped */
static int
mmCmpExch(medacc_state_t * state, const int i, const int j)
{
  return (mmless(state, i, j) && mmexchange(state , i, j));
}
 
/* maintains minheap property for all items below i/2. */
static void
minSortDown(medacc_state_t * state, int i)
{
  for (; i <= minCt(state); i *= 2)
    {
      if (i > 1 && i < minCt(state) && mmless(state, i + 1, i))
        ++i;

      if (!mmCmpExch(state, i, i / 2))
        break;
   }
}
 
/* maintains maxheap property for all items below i/2. (negative indexes) */
static void
maxSortDown(medacc_state_t * state, int i)
{
  for (; i >= -maxCt(state); i *= 2)
    {
      if (i < -1 && i > -maxCt(state) && mmless(state, i, i - 1))
        --i;

      if (!mmCmpExch(state, i / 2, i))
        break;
   }
}
 
/* maintains minheap property for all items above i, including median
   returns true if median changed */
static int
minSortUp(medacc_state_t * state, int i)
{
  while (i > 0 && mmCmpExch(state, i, i / 2))
    i /= 2;

  return (i == 0);
}
 
/* maintains maxheap property for all items above i, including median
   returns true if median changed */
static int
maxSortUp(medacc_state_t * state, int i)
{
  while (i<0 && mmCmpExch(state, i / 2, i))
    i /= 2;

  return (i == 0);
}

static const gsl_movstat_accum median_accum_type =
{
  medacc_size,
  medacc_init,
  medacc_insert,
  NULL, /* XXX FIXME */
  medacc_get
};

const gsl_movstat_accum *gsl_movstat_accum_median = &median_accum_type;
