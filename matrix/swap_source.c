/* matrix/swap_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
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
FUNCTION (gsl_matrix, swap_rows) (TYPE (gsl_matrix) * m,
                                 const size_t i, const size_t j)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (i >= size1)
    {
      GSL_ERROR ("first row index is out of range", GSL_EINVAL);
    }

  if (j >= size1)
    {
      GSL_ERROR ("second row index is out of range", GSL_EINVAL);
    }

  if (i != j)
    {
      ATOMIC *row1 = m->data + MULTIPLICITY * i * m->tda;
      ATOMIC *row2 = m->data + MULTIPLICITY * j * m->tda;
      
      size_t k;
      
      for (k = 0; k < MULTIPLICITY * size2; k++)
        {
          ATOMIC tmp = row1[k] ;
          row1[k] = row2[k] ;
          row2[k] = tmp ;
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_matrix, swap_columns) (TYPE (gsl_matrix) * m,
                                     const size_t i, const size_t j)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (i >= size2)
    {
      GSL_ERROR ("first column index is out of range", GSL_EINVAL);
    }

  if (j >= size2)
    {
      GSL_ERROR ("second column index is out of range", GSL_EINVAL);
    }

  if (i != j)
    {
      ATOMIC *col1 = m->data + MULTIPLICITY * i;
      ATOMIC *col2 = m->data + MULTIPLICITY * j;
      
      size_t p;
      
      for (p = 0; p < size1; p++)
        {
          size_t k;
          size_t n = p * MULTIPLICITY * m->tda;
 
          for (k = 0; k < MULTIPLICITY; k++)
            {
              ATOMIC tmp = col1[n+k] ;
              col1[n+k] = col2[n+k] ;
              col2[n+k] = tmp ;
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_matrix, swap_rowcol) (TYPE (gsl_matrix) * m,
                                    const size_t i, const size_t j)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (size1 != size2)
    {
      GSL_ERROR ("matrix must be square to swap row and column", GSL_ENOTSQR);
    }

  if (i >= size1)
    {
      GSL_ERROR ("row index is out of range", GSL_EINVAL);
    }

  if (j >= size2)
    {
      GSL_ERROR ("column index is out of range", GSL_EINVAL);
    }

  {
    ATOMIC *row = m->data + MULTIPLICITY * i * m->tda;
    ATOMIC *col = m->data + MULTIPLICITY * j;
      
    size_t p;
    
    for (p = 0; p < size1; p++)
      {
        size_t k;

        size_t r = p * MULTIPLICITY;
        size_t c = p * MULTIPLICITY * m->tda;
        
          for (k = 0; k < MULTIPLICITY; k++)
            {
              ATOMIC tmp = col[c+k] ;
              col[c+k] = row[r+k] ;
              row[r+k] = tmp ;
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_matrix, transpose) (TYPE (gsl_matrix) * m)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;
  size_t i, j, k;

  if (size1 != size2)
    {
      GSL_ERROR ("matrix must be square to take transpose", GSL_ENOTSQR);
    }

  for (i = 0; i < size1; i++)
    {
      for (j = i + 1 ; j < size2 ; j++) 
        {
          for (k = 0; k < MULTIPLICITY; k++)
            {
              size_t e1 = (i *  m->tda + j) * MULTIPLICITY + k ;
              size_t e2 = (j *  m->tda + i) * MULTIPLICITY + k ;
              {
                ATOMIC tmp = m->data[e1] ;
                m->data[e1] = m->data[e2] ;
                m->data[e2] = tmp ;
              }
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_matrix, transpose_memcpy) (TYPE (gsl_matrix) * dest, 
                                         const TYPE (gsl_matrix) * src)
{
  const size_t src_size1 = src->size1;
  const size_t src_size2 = src->size2;
  const size_t dest_size1 = dest->size1;
  const size_t dest_size2 = dest->size2;
  size_t i;

  if (dest_size2 != src_size1 || dest_size1 != src_size2)
    {
      GSL_ERROR ("dimensions of dest matrix must be transpose of src matrix", 
                 GSL_EBADLEN);
    }

#if defined(BASE_DOUBLE) || defined(BASE_FLOAT) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

  for (i = 0; i < src_size1; ++i)
    {
      VIEW (gsl_vector, const_view) a = FUNCTION (gsl_matrix, const_row) (src, i);
      VIEW (gsl_vector, view) b = FUNCTION (gsl_matrix, column) (dest, i);

#if defined(BASE_DOUBLE)
      gsl_blas_dcopy(&a.vector, &b.vector);
#elif defined(BASE_FLOAT)
      gsl_blas_scopy(&a.vector, &b.vector);
#elif defined(BASE_GSL_COMPLEX)
      gsl_blas_zcopy(&a.vector, &b.vector);
#elif defined(BASE_GSL_COMPLEX_FLOAT)
      gsl_blas_ccopy(&a.vector, &b.vector);
#endif
    }

#else

  for (i = 0; i < dest_size1; i++)
    {
      size_t j, k;

      for (j = 0 ; j < dest_size2; j++) 
        {
          for (k = 0; k < MULTIPLICITY; k++)
            {
              size_t e1 = (i *  dest->tda + j) * MULTIPLICITY + k ;
              size_t e2 = (j *  src->tda + i) * MULTIPLICITY + k ;

              dest->data[e1] = src->data[e2] ;
            }
        }
    }

#endif

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_matrix, transpose_tricpy) (CBLAS_UPLO_t Uplo_src,
                                         CBLAS_DIAG_t Diag,
                                         TYPE (gsl_matrix) * dest,
                                         const TYPE (gsl_matrix) * src)
{
  const size_t M = src->size1;
  const size_t N = src->size2;
  const size_t K = GSL_MIN(M, N);
  size_t i;

  if (M != dest->size2 || N != dest->size1)
    {
      GSL_ERROR ("matrix sizes are different", GSL_EBADLEN);
    }

#if defined(BASE_DOUBLE) || defined(BASE_FLOAT) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

  if (Uplo_src == CblasLower)
    {
      for (i = 1; i < K; i++)
        {
          VIEW (gsl_vector, const_view) a = FUNCTION (gsl_matrix, const_subrow) (src, i, 0, i);
          VIEW (gsl_vector, view) b = FUNCTION (gsl_matrix, subcolumn) (dest, i, 0, i);

#if defined(BASE_DOUBLE)
          gsl_blas_dcopy(&a.vector, &b.vector);
#elif defined(BASE_FLOAT)
          gsl_blas_scopy(&a.vector, &b.vector);
#elif defined(BASE_GSL_COMPLEX)
          gsl_blas_zcopy(&a.vector, &b.vector);
#elif defined(BASE_GSL_COMPLEX_FLOAT)
          gsl_blas_ccopy(&a.vector, &b.vector);
#endif
        }
    }
  else if (Uplo_src == CblasUpper)
    {
      for (i = 0; i < K - 1; i++)
        {
          VIEW (gsl_vector, const_view) a = FUNCTION (gsl_matrix, const_subrow) (src, i, i + 1, K - i - 1);
          VIEW (gsl_vector, view) b = FUNCTION (gsl_matrix, subcolumn) (dest, i, i + 1, K - i - 1);

#if defined(BASE_DOUBLE)
          gsl_blas_dcopy(&a.vector, &b.vector);
#elif defined(BASE_FLOAT)
          gsl_blas_scopy(&a.vector, &b.vector);
#elif defined(BASE_GSL_COMPLEX)
          gsl_blas_zcopy(&a.vector, &b.vector);
#elif defined(BASE_GSL_COMPLEX_FLOAT)
          gsl_blas_ccopy(&a.vector, &b.vector);
#endif
        }
    }

  if (Diag == CblasNonUnit)
    {
      VIEW (gsl_vector, const_view) a = FUNCTION (gsl_matrix, const_diagonal) (src);
      VIEW (gsl_vector, view) b = FUNCTION (gsl_matrix, diagonal) (dest);

#if defined(BASE_DOUBLE)
      gsl_blas_dcopy(&a.vector, &b.vector);
#elif defined(BASE_FLOAT)
      gsl_blas_scopy(&a.vector, &b.vector);
#elif defined(BASE_GSL_COMPLEX)
      gsl_blas_zcopy(&a.vector, &b.vector);
#elif defined(BASE_GSL_COMPLEX_FLOAT)
      gsl_blas_ccopy(&a.vector, &b.vector);
#endif
    }

#else

  {
    const size_t src_tda = src->tda ;
    const size_t dest_tda = dest->tda ;
    size_t j, k;

    if (Uplo_src == CblasLower)
      {
        /* copy lower triangle of src to upper triangle of dest */
        for (i = 0; i < K; i++)
          {
            for (j = 0; j < i; j++)
              {
                for (k = 0; k < MULTIPLICITY; k++)
                  {
                    size_t e1 = (j *  dest_tda + i) * MULTIPLICITY + k ;
                    size_t e2 = (i *  src_tda + j) * MULTIPLICITY + k ;
                    dest->data[e1] = src->data[e2];
                  }
              }
          }
      }
    else if (Uplo_src == CblasUpper)
      {
        /* copy upper triangle of src to lower triangle of dest */
        for (i = 0; i < K; i++)
          {
            for (j = i + 1; j < K; j++)
              {
                for (k = 0; k < MULTIPLICITY; k++)
                  {
                    size_t e1 = (j *  dest_tda + i) * MULTIPLICITY + k ;
                    size_t e2 = (i *  src_tda + j) * MULTIPLICITY + k ;
                    dest->data[e1] = src->data[e2];
                  }
              }
          }
      }
    else
      {
        GSL_ERROR ("invalid Uplo_src parameter", GSL_EINVAL);
      }

    if (Diag == CblasNonUnit)
      {
        for (i = 0; i < K; i++)
          {
            for (k = 0; k < MULTIPLICITY; k++)
              {
                size_t e1 = (i * dest_tda + i) * MULTIPLICITY + k ;
                size_t e2 = (i * src_tda + i) * MULTIPLICITY + k ;
                dest->data[e1] = src->data[e2];
              }
          }
      }
  }

#endif

  return GSL_SUCCESS;
}
