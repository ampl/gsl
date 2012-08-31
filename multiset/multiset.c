/* multiset/multiset.c
 * based on combination/combination.c by Szymon Jaroszewicz
 * based on permutation/permutation.c by Brian Gough
 *
 * Copyright (C) 2001 Szymon Jaroszewicz
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiset.h>

size_t
gsl_multiset_n (const gsl_multiset * c)
{
  return c->n ;
}

size_t
gsl_multiset_k (const gsl_multiset * c)
{
  return c->k ;
}

size_t *
gsl_multiset_data (const gsl_multiset * c)
{
  return c->data ;
}

int
gsl_multiset_valid (gsl_multiset * c)
{
  const size_t n = c->n ;
  const size_t k = c->k ;

  size_t i, j ;

  if( k > n )
    {
      GSL_ERROR("multiset has k greater than n", GSL_FAILURE) ;
    }
  for (i = 0; i < k; i++)
    {
      const size_t ci = c->data[i];

      if (ci >= n)
        {
          GSL_ERROR("multiset index outside range", GSL_FAILURE) ;
        }

      for (j = 0; j < i; j++)
        {
          if (c->data[j] > ci)
            {
              GSL_ERROR("multiset indices not in increasing order",
                        GSL_FAILURE) ;
            }
        }
    }

  return GSL_SUCCESS;
}


int
gsl_multiset_next (gsl_multiset * c)
{
  /* Replaces c with the next multiset (in the standard lexicographical
   * ordering).  Returns GSL_FAILURE if there is no next multiset.
   */
  const size_t n = c->n;
  const size_t k = c->k;
  size_t *data = c->data;
  size_t i;

  if(k == 0)
    {
      return GSL_FAILURE;
    }
  i = k - 1;

  while(i > 0 && data[i] == n-1)
    {
      --i;
    }

  if (i == 0 && data[0] == n-1)
    {
      return GSL_FAILURE;
    }

  ++data[i];

  while(i < k-1)
    {
      data[i+1] = data[i];
      ++i;
    }

  return GSL_SUCCESS;
}

int
gsl_multiset_prev (gsl_multiset * c)
{
  /* Replaces c with the previous multiset (in the standard
   * lexicographical ordering).  Returns GSL_FAILURE if there is no
   * previous multiset.
   */
  const size_t n = c->n;
  const size_t k = c->k;
  size_t *data = c->data;
  size_t i;

  if(k == 0)
    {
      return GSL_FAILURE;
    }
  i = k - 1;

  while(i > 0 && data[i-1] == data[i])
    {
      --i;
    }

  if(i == 0 && data[i] == 0)
    {
      return GSL_FAILURE;
    }

  data[i]--;

  if (data[i] < n-1)
    {
      while (i < k-1) {
        data[++i] = n - 1;
      }
    }

  return GSL_SUCCESS;
}

int
gsl_multiset_memcpy (gsl_multiset * dest, const gsl_multiset * src)
{
   const size_t src_n = src->n;
   const size_t src_k = src->k;
   const size_t dest_n = dest->n;
   const size_t dest_k = dest->k;

   if (src_n != dest_n || src_k != dest_k)
     {
       GSL_ERROR ("multiset lengths are not equal", GSL_EBADLEN);
     }

   {
     size_t j;

     for (j = 0; j < src_k; j++)
       {
         dest->data[j] = src->data[j];
       }
   }

   return GSL_SUCCESS;
}
