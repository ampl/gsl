/* spmatrix/copy_source.c
 * 
 * Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Patrick Alken
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

int
FUNCTION (gsl_spmatrix, memcpy) (TYPE (gsl_spmatrix) * dest, const TYPE (gsl_spmatrix) * src)
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
      int status = GSL_SUCCESS;
      size_t n, r;

      if (dest->nzmax < src->nz)
        {
          status = FUNCTION (gsl_spmatrix, realloc) (src->nz, dest);
          if (status)
            return status;
        }

      /* copy src arrays to dest */
      if (GSL_SPMATRIX_ISCOO(src))
        {
          void *ptr;

          for (n = 0; n < src->nz; ++n)
            {
              dest->i[n] = src->i[n];
              dest->p[n] = src->p[n];

              for (r = 0; r < MULTIPLICITY; ++r)
                dest->data[MULTIPLICITY * n + r] = src->data[MULTIPLICITY * n + r];

              /* copy binary tree data */
              ptr = gsl_bst_insert(&dest->data[MULTIPLICITY * n], dest->tree);
              if (ptr != NULL)
                {
                  GSL_ERROR("detected duplicate entry", GSL_EINVAL);
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSC(src))
        {
          for (n = 0; n < src->nz; ++n)
            {
              dest->i[n] = src->i[n];

              for (r = 0; r < MULTIPLICITY; ++r)
                dest->data[MULTIPLICITY * n + r] = src->data[MULTIPLICITY * n + r];
            }

          for (n = 0; n < src->size2 + 1; ++n)
            {
              dest->p[n] = src->p[n];
            }
        }
      else if (GSL_SPMATRIX_ISCSR(src))
        {
          for (n = 0; n < src->nz; ++n)
            {
              dest->i[n] = src->i[n];

              for (r = 0; r < MULTIPLICITY; ++r)
                dest->data[MULTIPLICITY * n + r] = src->data[MULTIPLICITY * n + r];
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

      return status;
    }
}
