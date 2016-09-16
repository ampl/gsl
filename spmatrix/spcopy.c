/* spcopy.c
 * 
 * Copyright (C) 2014 Patrick Alken
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
#include <math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_errno.h>

#include "avl.c"

int
gsl_spmatrix_memcpy(gsl_spmatrix *dest, const gsl_spmatrix *src)
{
  const size_t M = src->size1;
  const size_t N = src->size2;

  if (M != dest->size1 || N != dest->size2)
    {
      GSL_ERROR("matrix sizes are different", GSL_EBADLEN);
    }
  else if (dest->sptype != src->sptype)
    {
      GSL_ERROR("cannot copy matrices of different storage formats",
                GSL_EINVAL);
    }
  else
    {
      int s = GSL_SUCCESS;
      size_t n;

      if (dest->nzmax < src->nz)
        {
          s = gsl_spmatrix_realloc(src->nz, dest);
          if (s)
            return s;
        }

      /* copy indices and data to dest */
      if (GSL_SPMATRIX_ISTRIPLET(src))
        {
          void *ptr;

          for (n = 0; n < src->nz; ++n)
            {
              dest->i[n] = src->i[n];
              dest->p[n] = src->p[n];
              dest->data[n] = src->data[n];

              /* copy binary tree data */
              ptr = avl_insert(dest->tree_data->tree, &dest->data[n]);
              if (ptr != NULL)
                {
                  GSL_ERROR("detected duplicate entry", GSL_EINVAL);
                }
            }
        }
      else if (GSL_SPMATRIX_ISCCS(src))
        {
          for (n = 0; n < src->nz; ++n)
            {
              dest->i[n] = src->i[n];
              dest->data[n] = src->data[n];
            }

          for (n = 0; n < src->size2 + 1; ++n)
            {
              dest->p[n] = src->p[n];
            }
        }
      else if (GSL_SPMATRIX_ISCRS(src))
        {
          for (n = 0; n < src->nz; ++n)
            {
              dest->i[n] = src->i[n];
              dest->data[n] = src->data[n];
            }

          for (n = 0; n < src->size1 + 1; ++n)
            {
              dest->p[n] = src->p[n];
            }
        }
      else
        {
          GSL_ERROR("invalid matrix type for src", GSL_EINVAL);
        }

      dest->nz = src->nz;

      return s;
    }
} /* gsl_spmatrix_memcpy() */
