/* movstat/ringbuf.c
 *
 * Ring buffer module
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

#ifndef __GSL_RINGBUF_C__
#define __GSL_RINGBUF_C__

/*typedef int ringbuf_type;*/

typedef struct
{
  ringbuf_type_t *array;
  int head;
  int tail;
  int size;       /* total elements allocated */
} ringbuf;

static size_t ringbuf_size(const size_t n);
static int ringbuf_empty(ringbuf * d);
static int ringbuf_is_empty(const ringbuf * d);
static int ringbuf_is_full(const ringbuf * d);
static int ringbuf_insert(const ringbuf_type_t x, ringbuf * d);
static int ringbuf_pop_back(ringbuf * b);
static ringbuf_type_t ringbuf_peek(const int i, const ringbuf * b);
static ringbuf_type_t ringbuf_peek_front(const ringbuf * d);
static ringbuf_type_t ringbuf_peek_back(const ringbuf * d);
static size_t ringbuf_copy(double * dest, const ringbuf * b);
static int ringbuf_n(const ringbuf * b);

static size_t
ringbuf_size(const size_t n)
{
  size_t size = 0;

  size += sizeof(ringbuf);
  size += n * sizeof(ringbuf_type_t); /* b->array */

  return size;
}

static int
ringbuf_init(const size_t n, ringbuf * b)
{
  b->array = (ringbuf_type_t *) ((char *) b + sizeof(ringbuf));
  b->head = -1;
  b->tail = 0;
  b->size = (int) n;
  return GSL_SUCCESS;
}

/* empty the buffer */
static int
ringbuf_empty(ringbuf * b)
{
  b->head = -1;
  b->tail = 0;
  return GSL_SUCCESS;
}

/* check if buffer is empty */
static int
ringbuf_is_empty(const ringbuf * b)
{
  return (b->head == -1);
}

/* check if buffer is full */
static int
ringbuf_is_full(const ringbuf * b)
{
  return ((b->head == 0 && b->tail == b->size - 1) ||
          (b->head == b->tail + 1));
}

/* insert element into buffer, overwriting oldest element if necessary */
static int
ringbuf_insert(const ringbuf_type_t x, ringbuf * b)
{
  if (b->head == -1)     /* buffer is empty */
    {
      b->head = 0;
      b->tail = 0;
    }
  else if (b->head == 0) /* head is in first position, wrap to end */
    {
      b->head = b->size - 1;

      if (b->tail == b->head && b->size > 1)
        --(b->tail);     /* buffer is full so decrease tail */
    }
  else                   /* decrement head */
    {
      --(b->head);

      if (b->tail == b->head)
        {
          /* buffer is full so update tail */
          if (b->tail == 0)
            b->tail = b->size - 1;
          else
            --(b->tail);
        }
    }

  /* insert element */
  b->array[b->head] = x;

  return GSL_SUCCESS;
}

static int
ringbuf_pop_back(ringbuf * b)
{
  if (ringbuf_is_empty(b) || b->tail < 0)
    {
      GSL_ERROR("buffer is empty", GSL_EBADLEN);
    }
  else
    {
      if (b->head == b->tail) /* buffer has only one element */
        {
          b->head = -1;
          b->tail = -1;
        }
      else if (b->tail == 0)  /* tail is in first position, wrap to end */
        {
          b->tail = b->size - 1;
        }
      else                    /* decrement tail */
        {
          --(b->tail);
        }

      return GSL_SUCCESS;
    }
}

static ringbuf_type_t
ringbuf_peek(const int i, const ringbuf * b)
{
  if (ringbuf_is_empty(b))
    {
      GSL_ERROR("buffer is empty", GSL_EBADLEN);
    }
  else
    {
      return b->array[(b->head + i) % b->size];
    }
}

static ringbuf_type_t
ringbuf_peek_front(const ringbuf * b)
{
  if (ringbuf_is_empty(b))
    {
      GSL_ERROR("buffer is empty", GSL_EBADLEN);
    }
  else
    {
      return b->array[b->head];
    }
}

static ringbuf_type_t
ringbuf_peek_back(const ringbuf * b)
{
  if (ringbuf_is_empty(b) || b->tail < 0)
    {
      GSL_ERROR("buffer is empty", GSL_EBADLEN);
    }
  else
    {
      return b->array[b->tail];
    }
}

static size_t
ringbuf_copy(double * dest, const ringbuf * b)
{
  if (ringbuf_is_empty(b) || b->tail < 0)
    {
      return 0;
    }
  else
    {
      const int n = ringbuf_n(b);
      int i;

      for (i = 0; i < n; ++i)
        dest[i] = b->array[(b->head + i) % b->size];

      return (size_t) n;
    }
}

static int
ringbuf_n(const ringbuf * b)
{
  const int n = (b->head > b->tail) ? (b->size - b->head + b->tail + 1) : (b->tail - b->head + 1);
  return n;
}

#endif /* __GSL_RINGBUF_C__ */
