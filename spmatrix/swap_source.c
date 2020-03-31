/* spmatrix/swap_source.c
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
FUNCTION (gsl_spmatrix, transpose) (TYPE (gsl_spmatrix) * m)
{
  /* swap dimensions - this must be done before the tree_rebuild step */
  if (m->size1 != m->size2)
    {
      size_t tmp = m->size1;
      m->size1 = m->size2;
      m->size2 = tmp;
    }
  
  if (GSL_SPMATRIX_ISCOO(m))
    {
      size_t n;

      /* swap row/column indices */
      for (n = 0; n < m->nz; ++n)
        {
          int tmp = m->p[n];
          m->p[n] = m->i[n];
          m->i[n] = tmp;
        }

      /* need to rebuild binary tree, or element searches won't
       * work correctly with transposed indices */
      FUNCTION (gsl_spmatrix, tree_rebuild) (m);
    }
  else if (GSL_SPMATRIX_ISCSC(m))
    {
      m->sptype = GSL_SPMATRIX_CSR;
    }
  else if (GSL_SPMATRIX_ISCSR(m))
    {
      m->sptype = GSL_SPMATRIX_CSC;
    }
  else
    {
      GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
    }
  
  return GSL_SUCCESS;
}

/* FIXME XXX: deprecated function */
int
FUNCTION (gsl_spmatrix, transpose2) (TYPE (gsl_spmatrix) * m)
{
  return FUNCTION (gsl_spmatrix, transpose) (m);
}

int
FUNCTION (gsl_spmatrix, transpose_memcpy) (TYPE (gsl_spmatrix) * dest, const TYPE (gsl_spmatrix) * src)
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
      int status = GSL_SUCCESS;
      const size_t nz = src->nz;

      if (dest->nzmax < src->nz)
        {
          status = FUNCTION (gsl_spmatrix, realloc) (src->nz, dest);
          if (status)
            return status;
        }

      if (GSL_SPMATRIX_ISCOO(src))
        {
          size_t n, r;
          void *ptr;

          for (n = 0; n < nz; ++n)
            {
              dest->i[n] = src->p[n];
              dest->p[n] = src->i[n];

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
          int * Ai = src->i;
          int * Ap = src->p;
          ATOMIC * Ad = src->data;
          int * ATi = dest->i;
          int * ATp = dest->p;
          ATOMIC * ATd = dest->data;
          int * w = dest->work.work_int;
          int p;
          size_t j, r;

          /* initialize to 0 */
          for (j = 0; j < M + 1; ++j)
            ATp[j] = 0;

          /* compute row counts of A (= column counts for A^T) */
          for (j = 0; j < nz; ++j)
            ATp[Ai[j]]++;

          /* compute row pointers for A (= column pointers for A^T) */
          gsl_spmatrix_cumsum(M, ATp);

          /* make copy of row pointers */
          for (j = 0; j < M; ++j)
            w[j] = ATp[j];

          for (j = 0; j < N; ++j)
            {
              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  int k = w[Ai[p]]++;
                  ATi[k] = j;

                  for (r = 0; r < MULTIPLICITY; ++r)
                    ATd[MULTIPLICITY * k + r] = Ad[MULTIPLICITY * p + r];
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSR(src))
        {
          int * Aj = src->i;
          int * Ap = src->p;
          ATOMIC * Ad = src->data;
          int * ATj = dest->i;
          int * ATp = dest->p;
          ATOMIC * ATd = dest->data;
          int * w = dest->work.work_int;
          int p;
          size_t i, r;

          /* initialize to 0 */
          for (i = 0; i < N + 1; ++i)
            ATp[i] = 0;

          /* compute column counts of A (= row counts for A^T) */
          for (i = 0; i < nz; ++i)
            ATp[Aj[i]]++;

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

                  for (r = 0; r < MULTIPLICITY; ++r)
                    ATd[MULTIPLICITY * k + r] = Ad[MULTIPLICITY * p + r];
                }
            }
        }
      else
        {
          GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
        }

      dest->nz = nz;

      return status;
    }
}
