/* spcompress.c
 * 
 * Copyright (C) 2012-2014, 2016 Patrick Alken
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>

/*
gsl_spmatrix_ccs()
  Create a sparse matrix in compressed column format

Inputs: T - sparse matrix in triplet format

Return: pointer to new matrix (should be freed when finished with it)
*/

gsl_spmatrix *
gsl_spmatrix_ccs(const gsl_spmatrix *T)
{
  if (!GSL_SPMATRIX_ISTRIPLET(T))
    {
      GSL_ERROR_NULL("matrix must be in triplet format", GSL_EINVAL);
    }
  else
    {
      const size_t *Tj; /* column indices of triplet matrix */
      size_t *Cp;       /* column pointers of compressed column matrix */
      size_t *w;        /* copy of column pointers */
      gsl_spmatrix *m;
      size_t n;

      m = gsl_spmatrix_alloc_nzmax(T->size1, T->size2, T->nz,
                                   GSL_SPMATRIX_CCS);
      if (!m)
        return NULL;

      Tj = T->p;
      Cp = m->p;

      /* initialize column pointers to 0 */
      for (n = 0; n < m->size2 + 1; ++n)
        Cp[n] = 0;

      /*
       * compute the number of elements in each column:
       * Cp[j] = # non-zero elements in column j
       */
      for (n = 0; n < T->nz; ++n)
        Cp[Tj[n]]++;

      /* compute column pointers: p[j] = p[j-1] + nnz[j-1] */
      gsl_spmatrix_cumsum(m->size2, Cp);

      /* make a copy of the column pointers */
      w = m->work_sze;
      for (n = 0; n < m->size2; ++n)
        w[n] = Cp[n];

      /* transfer data from triplet format to CCS */
      for (n = 0; n < T->nz; ++n)
        {
          size_t k = w[Tj[n]]++;
          m->i[k] = T->i[n];
          m->data[k] = T->data[n];
        }

      m->nz = T->nz;

      return m;
    }
}

gsl_spmatrix *
gsl_spmatrix_compcol(const gsl_spmatrix *T)
{
  return gsl_spmatrix_ccs(T);
}

/*
gsl_spmatrix_crs()
  Create a sparse matrix in compressed row format

Inputs: T - sparse matrix in triplet format

Return: pointer to new matrix (should be freed when finished with it)
*/

gsl_spmatrix *
gsl_spmatrix_crs(const gsl_spmatrix *T)
{
  if (!GSL_SPMATRIX_ISTRIPLET(T))
    {
      GSL_ERROR_NULL("matrix must be in triplet format", GSL_EINVAL);
    }
  else
    {
      const size_t *Ti; /* row indices of triplet matrix */
      size_t *Cp;       /* row pointers of compressed row matrix */
      size_t *w;        /* copy of column pointers */
      gsl_spmatrix *m;
      size_t n;

      m = gsl_spmatrix_alloc_nzmax(T->size1, T->size2, T->nz,
                                   GSL_SPMATRIX_CRS);
      if (!m)
        return NULL;

      Ti = T->i;
      Cp = m->p;

      /* initialize row pointers to 0 */
      for (n = 0; n < m->size1 + 1; ++n)
        Cp[n] = 0;

      /*
       * compute the number of elements in each row:
       * Cp[i] = # non-zero elements in row i
       */
      for (n = 0; n < T->nz; ++n)
        Cp[Ti[n]]++;

      /* compute row pointers: p[i] = p[i-1] + nnz[i-1] */
      gsl_spmatrix_cumsum(m->size1, Cp);

      /* make a copy of the row pointers */
      w = m->work_sze;
      for (n = 0; n < m->size1; ++n)
        w[n] = Cp[n];

      /* transfer data from triplet format to CRS */
      for (n = 0; n < T->nz; ++n)
        {
          size_t k = w[Ti[n]]++;
          m->i[k] = T->p[n];
          m->data[k] = T->data[n];
        }

      m->nz = T->nz;

      return m;
    }
}

/*
gsl_spmatrix_cumsum()

Compute the cumulative sum:

p[j] = Sum_{k=0...j-1} c[k]

0 <= j < n + 1

Alternatively,
p[0] = 0
p[j] = p[j - 1] + c[j - 1]

Inputs: n - length of input array
        c - (input/output) array of size n + 1
            on input, contains the n values c[k]
            on output, contains the n + 1 values p[j]

Return: success or error
*/

void
gsl_spmatrix_cumsum(const size_t n, size_t *c)
{
  size_t sum = 0;
  size_t k;

  for (k = 0; k < n; ++k)
    {
      size_t ck = c[k];
      c[k] = sum;
      sum += ck;
    }

  c[n] = sum;
} /* gsl_spmatrix_cumsum() */
