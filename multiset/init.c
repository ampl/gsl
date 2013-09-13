/* multiset/init.c
 * based on combination/init.c by Szymon Jaroszewicz
 * based on permutation/init.c by Brian Gough
 *
 * Copyright (C) 2001 Szymon Jaroszewicz
 * Copyright (C) 2009 Brian Gough
 * Copyright (C) 2009 Rhys Ulerich
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
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiset.h>

gsl_multiset *
gsl_multiset_alloc (const size_t n, const size_t k)
{
  gsl_multiset * c;

  if (n == 0)
    {
      GSL_ERROR_VAL ("multiset parameter n must be positive integer",
                        GSL_EDOM, 0);
    }
  c = (gsl_multiset *) malloc (sizeof (gsl_multiset));

  if (c == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for multiset struct",
                        GSL_ENOMEM, 0);
    }

  if (k > 0)
    {
      c->data = (size_t *) malloc (k * sizeof (size_t));

      if (c->data == 0)
        {
          free (c);             /* exception in constructor, avoid memory leak */

          GSL_ERROR_VAL ("failed to allocate space for multiset data",
                         GSL_ENOMEM, 0);
        }
    }
  else
    {
      c->data = 0;
    }

  c->n = n;
  c->k = k;

  return c;
}

gsl_multiset *
gsl_multiset_calloc (const size_t n, const size_t k)
{
  size_t i;

  gsl_multiset * c =  gsl_multiset_alloc (n, k);

  if (c == 0)
    return 0;

  /* initialize multiset to repeated first element */

  for (i = 0; i < k; i++)
    {
      c->data[i] = 0;
    }

  return c;
}

void
gsl_multiset_init_first (gsl_multiset * c)
{
  const size_t k = c->k ;
  size_t i;

  /* initialize multiset to repeated first element */

  for (i = 0; i < k; i++)
    {
      c->data[i] = 0;
    }
}

void
gsl_multiset_init_last (gsl_multiset * c)
{
  const size_t k = c->k ;
  size_t i;
  size_t n = c->n;

  /* initialize multiset to repeated last element */

  for (i = 0; i < k; i++)
    {
      c->data[i] = n - 1;
    }
}

void
gsl_multiset_free (gsl_multiset * c)
{
  RETURN_IF_NULL (c);
  if (c->k > 0) free (c->data);
  free (c);
}
