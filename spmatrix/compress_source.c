/* spmatrix/compress_source.c
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

/*
gsl_spmatrix_csc()
  Create a sparse matrix in compressed column format

Inputs: T - sparse matrix in triplet format

Return: pointer to new matrix (should be freed when finished with it)
*/

int
FUNCTION (gsl_spmatrix, csc) (TYPE (gsl_spmatrix) * dest, const TYPE (gsl_spmatrix) * src)
{
  if (!GSL_SPMATRIX_ISCOO(src))
    {
      GSL_ERROR_NULL("input matrix must be in COO format", GSL_EINVAL);
    }
  else if (!GSL_SPMATRIX_ISCSC(dest))
    {
      GSL_ERROR_NULL("output matrix must be in CSC format", GSL_EINVAL);
    }
  else if (src->size1 != dest->size1 || src->size2 != dest->size2)
    {
      GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
    }
  else
    {
      int status;
      const int *Tj = src->p; /* column indices of triplet matrix */
      int *Cp;                /* column pointers of compressed column matrix */
      int *w;                 /* copy of column pointers */
      size_t n, r;

      if (dest->nzmax < src->nz)
        {
          status = FUNCTION (gsl_spmatrix, realloc) (src->nz, dest);
          if (status)
            return status;
        }

      Cp = dest->p;

      /* initialize column pointers to 0 */
      for (n = 0; n < dest->size2 + 1; ++n)
        Cp[n] = 0;

      /*
       * compute the number of elements in each column:
       * Cp[j] = # non-zero elements in column j
       */
      for (n = 0; n < src->nz; ++n)
        Cp[Tj[n]]++;

      /* compute column pointers: p[j] = p[j-1] + nnz[j-1] */
      gsl_spmatrix_cumsum(dest->size2, Cp);

      /* make a copy of the column pointers */
      w = dest->work.work_int;
      for (n = 0; n < dest->size2; ++n)
        w[n] = Cp[n];

      /* transfer data from triplet format to CSC */
      for (n = 0; n < src->nz; ++n)
        {
          int k = w[Tj[n]]++;
          dest->i[k] = src->i[n];

          for (r = 0; r < MULTIPLICITY; ++r)
            dest->data[MULTIPLICITY * k + r] = src->data[MULTIPLICITY * n + r];
        }

      dest->nz = src->nz;

      return GSL_SUCCESS;
    }
}

/* XXX deprecated function */
TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, compcol) (const TYPE (gsl_spmatrix) * src)
{
  return FUNCTION (gsl_spmatrix, ccs) (src);
}

/* XXX deprecated function */
TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, ccs) (const TYPE (gsl_spmatrix) * src)
{
  TYPE (gsl_spmatrix) * dest = FUNCTION (gsl_spmatrix, alloc_nzmax) (src->size1, src->size2, src->nz, GSL_SPMATRIX_CSC);
  FUNCTION (gsl_spmatrix, csc) (dest, src);
  return dest;
}

/*
gsl_spmatrix_csr()
  Create a sparse matrix in compressed row format

Inputs: dest - (output) sparse matrix in CSR format
        src  - sparse matrix in triplet format

Return: success/error
*/

int
FUNCTION (gsl_spmatrix, csr) (TYPE (gsl_spmatrix) * dest, const TYPE (gsl_spmatrix) * src)
{
  if (!GSL_SPMATRIX_ISCOO(src))
    {
      GSL_ERROR("input matrix must be in COO format", GSL_EINVAL);
    }
  else if (!GSL_SPMATRIX_ISCSR(dest))
    {
      GSL_ERROR("output matrix must be in CSR format", GSL_EINVAL);
    }
  else if (src->size1 != dest->size1 || src->size2 != dest->size2)
    {
      GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
    }
  else
    {
      int status;
      const int *Ti = src->i; /* row indices of triplet matrix */
      int *Cp;                /* row pointers of compressed row matrix */
      int *w;                 /* copy of column pointers */
      size_t n, r;

      if (dest->nzmax < src->nz)
        {
          status = FUNCTION (gsl_spmatrix, realloc) (src->nz, dest);
          if (status)
            return status;
        }

      Cp = dest->p;

      /* initialize row pointers to 0 */
      for (n = 0; n < dest->size1 + 1; ++n)
        Cp[n] = 0;

      /*
       * compute the number of elements in each row:
       * Cp[i] = # non-zero elements in row i
       */
      for (n = 0; n < src->nz; ++n)
        Cp[Ti[n]]++;

      /* compute row pointers: p[i] = p[i-1] + nnz[i-1] */
      gsl_spmatrix_cumsum(dest->size1, Cp);

      /* make a copy of the row pointers */
      w = dest->work.work_int;
      for (n = 0; n < dest->size1; ++n)
        w[n] = Cp[n];

      /* transfer data from triplet format to CSR */
      for (n = 0; n < src->nz; ++n)
        {
          int k = w[Ti[n]]++;
          dest->i[k] = src->p[n];

          for (r = 0; r < MULTIPLICITY; ++r)
            dest->data[MULTIPLICITY * k + r] = src->data[MULTIPLICITY * n + r];
        }

      dest->nz = src->nz;

      return GSL_SUCCESS;
    }
}

/* XXX deprecated function */
TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, crs) (const TYPE (gsl_spmatrix) * src)
{
  TYPE (gsl_spmatrix) * dest = FUNCTION (gsl_spmatrix, alloc_nzmax) (src->size1, src->size2, src->nz, GSL_SPMATRIX_CSR);
  FUNCTION (gsl_spmatrix, csr) (dest, src);
  return dest;
}

TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, compress) (const TYPE (gsl_spmatrix) * src, const int sptype)
{
  int status = GSL_SUCCESS;
  TYPE (gsl_spmatrix) * dest = FUNCTION (gsl_spmatrix, alloc_nzmax) (src->size1, src->size2, src->nz, sptype);

  if (dest == NULL)
    return NULL;

  if (sptype == GSL_SPMATRIX_CSC)
    {
      status = FUNCTION (gsl_spmatrix, csc) (dest, src);
    }
  else if (sptype == GSL_SPMATRIX_CSR)
    {
      status = FUNCTION (gsl_spmatrix, csr) (dest, src);
    }
  else if (sptype == GSL_SPMATRIX_COO)
    {
      /* make a copy of src */
      status = FUNCTION (gsl_spmatrix, memcpy) (dest, src);
    }
  else
    {
      status = GSL_EINVAL;
      GSL_ERROR_NULL ("unknown sparse matrix format", status);
    }

  if (status != GSL_SUCCESS)
    {
      FUNCTION (gsl_spmatrix, free) (dest);
      return NULL;
    }

  return dest;
}
