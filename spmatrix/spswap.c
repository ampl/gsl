/* spswap.c
 * 
 * Copyright (C) 2014 Patrick Alken
 * Copyright (C) 2016 Alexis Tantet
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>

#include "avl.c"

/*
gsl_spmatrix_transpose()
  Replace the sparse matrix src by its transpose,
keeping the matrix in the same storage format

Inputs: A - (input/output) sparse matrix to transpose.
*/

int
gsl_spmatrix_transpose(gsl_spmatrix * m)
{
  if (GSL_SPMATRIX_ISTRIPLET(m))
    {
      size_t n;

      /* swap row/column indices */
      for (n = 0; n < m->nz; ++n)
        {
          size_t tmp = m->p[n];
          m->p[n] = m->i[n];
          m->i[n] = tmp;
        }

      /* need to rebuild AVL tree, or element searches won't
       * work correctly with transposed indices */
      gsl_spmatrix_tree_rebuild(m);
    }
  else
    {
      GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
    }
  
  /* swap dimensions */
  if (m->size1 != m->size2)
    {
      size_t tmp = m->size1;
      m->size1 = m->size2;
      m->size2 = tmp;
    }
  
  return GSL_SUCCESS;
}

/*
gsl_spmatrix_transpose2()
  Replace the sparse matrix src by its transpose either by
  swapping its row and column indices if it is in triplet storage,
  or by switching its major if it is in compressed storage.

Inputs: A - (input/output) sparse matrix to transpose.
*/

int
gsl_spmatrix_transpose2(gsl_spmatrix * m)
{
  if (GSL_SPMATRIX_ISTRIPLET(m))
    {
      return gsl_spmatrix_transpose(m);
    }
  else if (GSL_SPMATRIX_ISCCS(m))
    {
      m->sptype = GSL_SPMATRIX_CRS;
    }
  else if (GSL_SPMATRIX_ISCRS(m))
    {
      m->sptype = GSL_SPMATRIX_CCS;
    }
  else
    {
      GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
    }
  
  /* swap dimensions */
  if (m->size1 != m->size2)
    {
      size_t tmp = m->size1;
      m->size1 = m->size2;
      m->size2 = tmp;
    }
  
  return GSL_SUCCESS;
}

int
gsl_spmatrix_transpose_memcpy(gsl_spmatrix *dest, const gsl_spmatrix *src)
{
  const size_t M = src->size1;
  const size_t N = src->size2;

  if (M != dest->size2 || N != dest->size1)
    {
      GSL_ERROR("dimensions of dest must be transpose of src matrix",
                GSL_EBADLEN);
    }
  else if (dest->sptype != src->sptype)
    {
      GSL_ERROR("cannot copy matrices of different storage formats",
                GSL_EINVAL);
    }
  else
    {
      int s = GSL_SUCCESS;
      const size_t nz = src->nz;

      if (dest->nzmax < src->nz)
        {
          s = gsl_spmatrix_realloc(src->nz, dest);
          if (s)
            return s;
        }

      if (GSL_SPMATRIX_ISTRIPLET(src))
        {
          size_t n;
          void *ptr;

          for (n = 0; n < nz; ++n)
            {
              dest->i[n] = src->p[n];
              dest->p[n] = src->i[n];
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
          size_t *Ai = src->i;
          size_t *Ap = src->p;
          double *Ad = src->data;
          size_t *ATi = dest->i;
          size_t *ATp = dest->p;
          double *ATd = dest->data;
          size_t *w = (size_t *) dest->work;
          size_t p, j;

          /* initialize to 0 */
          for (p = 0; p < M + 1; ++p)
            ATp[p] = 0;

          /* compute row counts of A (= column counts for A^T) */
          for (p = 0; p < nz; ++p)
            ATp[Ai[p]]++;

          /* compute row pointers for A (= column pointers for A^T) */
          gsl_spmatrix_cumsum(M, ATp);

          /* make copy of row pointers */
          for (j = 0; j < M; ++j)
            w[j] = ATp[j];

          for (j = 0; j < N; ++j)
            {
              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  size_t k = w[Ai[p]]++;
                  ATi[k] = j;
                  ATd[k] = Ad[p];
                }
            }
        }
      else if (GSL_SPMATRIX_ISCRS(src))
        {
          size_t *Aj = src->i;
          size_t *Ap = src->p;
          double *Ad = src->data;
          size_t *ATj = dest->i;
          size_t *ATp = dest->p;
          double *ATd = dest->data;
          size_t *w = (size_t *) dest->work;
          size_t p, i;

          /* initialize to 0 */
          for (p = 0; p < N + 1; ++p)
            ATp[p] = 0;

          /* compute column counts of A (= row counts for A^T) */
          for (p = 0; p < nz; ++p)
            ATp[Aj[p]]++;

          /* compute column pointers for A (= row pointers for A^T) */
          gsl_spmatrix_cumsum(N, ATp);

          /* make copy of column pointers */
          for (i = 0; i < N; ++i)
            w[i] = ATp[i];

          for (i = 0; i < M; ++i)
            {
              for (p = Ap[i]; p < Ap[i + 1]; ++p)
                {
                  size_t k = w[Aj[p]]++;
                  ATj[k] = i;
                  ATd[k] = Ad[p];
                }
            }
        }
      else
        {
          GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
        }

      dest->nz = nz;

      return s;
    }
} /* gsl_spmatrix_transpose_memcpy() */
