/* matrix/copy_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Copyright (C) 2019 Patrick Alken
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
FUNCTION (gsl_matrix, memcpy) (TYPE (gsl_matrix) * dest,
                               const TYPE (gsl_matrix) * src)
{
  const size_t src_size1 = src->size1;
  const size_t src_size2 = src->size2;
  const size_t dest_size1 = dest->size1;
  const size_t dest_size2 = dest->size2;
  size_t i;

  if (src_size1 != dest_size1 || src_size2 != dest_size2)
    {
      GSL_ERROR ("matrix sizes are different", GSL_EBADLEN);
    }

#if defined(BASE_DOUBLE) || defined(BASE_FLOAT) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

  for (i = 0; i < src_size1; ++i)
    {
      VIEW (gsl_vector, const_view) sv = FUNCTION (gsl_matrix, const_row) (src, i);
      VIEW (gsl_vector, view) dv = FUNCTION (gsl_matrix, row) (dest, i);

#if defined(BASE_DOUBLE)
      gsl_blas_dcopy(&sv.vector, &dv.vector);
#elif defined(BASE_FLOAT)
      gsl_blas_scopy(&sv.vector, &dv.vector);
#elif defined(BASE_GSL_COMPLEX)
      gsl_blas_zcopy(&sv.vector, &dv.vector);
#elif defined(BASE_GSL_COMPLEX_FLOAT)
      gsl_blas_ccopy(&sv.vector, &dv.vector);
#endif
    }

#else

  {
    const size_t src_tda = src->tda ;
    const size_t dest_tda = dest->tda ;
    size_t j;

    for (i = 0; i < src_size1 ; i++)
      {
        for (j = 0; j < MULTIPLICITY * src_size2; j++)
          {
            dest->data[MULTIPLICITY * dest_tda * i + j] 
              = src->data[MULTIPLICITY * src_tda * i + j];
          }
      }
  }

#endif

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_matrix, swap) (TYPE (gsl_matrix) * dest, TYPE (gsl_matrix) * src)
{
  const size_t src_size1 = src->size1;
  const size_t src_size2 = src->size2;
  const size_t dest_size1 = dest->size1;
  const size_t dest_size2 = dest->size2;
  size_t i;

  if (src_size1 != dest_size1 || src_size2 != dest_size2)
    {
      GSL_ERROR ("matrix sizes are different", GSL_EBADLEN);
    }

#if defined(BASE_DOUBLE) || defined(BASE_FLOAT) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

  for (i = 0; i < src_size1; ++i)
    {
      VIEW (gsl_vector, view) sv = FUNCTION (gsl_matrix, row) (src, i);
      VIEW (gsl_vector, view) dv = FUNCTION (gsl_matrix, row) (dest, i);

#if defined(BASE_DOUBLE)
      gsl_blas_dswap(&sv.vector, &dv.vector);
#elif defined(BASE_FLOAT)
      gsl_blas_sswap(&sv.vector, &dv.vector);
#elif defined(BASE_GSL_COMPLEX)
      gsl_blas_zswap(&sv.vector, &dv.vector);
#elif defined(BASE_GSL_COMPLEX_FLOAT)
      gsl_blas_cswap(&sv.vector, &dv.vector);
#endif
    }

#else

  {
    const size_t src_tda = src->tda ;
    const size_t dest_tda = dest->tda ;
    size_t j;

    for (i = 0; i < src_size1 ; i++)
      {
        for (j = 0; j < MULTIPLICITY * src_size2; j++)
          {
            ATOMIC tmp = src->data[MULTIPLICITY * src_tda * i + j];
            src->data[MULTIPLICITY * src_tda * i + j] 
              = dest->data[MULTIPLICITY * dest_tda * i + j];
            dest->data[MULTIPLICITY * dest_tda * i + j] = tmp ;
          }
      }
  }

#endif

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_matrix, tricpy) (CBLAS_UPLO_t Uplo,
                               CBLAS_DIAG_t Diag,
                               TYPE (gsl_matrix) * dest,
                               const TYPE (gsl_matrix) * src)
{
  const size_t M = src->size1;
  const size_t N = src->size2;
  size_t i;

  if (M != dest->size1 || N != dest->size2)
    {
      GSL_ERROR ("matrix sizes are different", GSL_EBADLEN);
    }

#if defined(BASE_DOUBLE) || defined(BASE_FLOAT) || defined(BASE_GSL_COMPLEX) || defined(BASE_GSL_COMPLEX_FLOAT)

  if (Uplo == CblasLower)
    {
      for (i = 1; i < M; i++)
        {
          size_t k = GSL_MIN (i, N);
          VIEW (gsl_vector, const_view) a = FUNCTION (gsl_matrix, const_subrow) (src, i, 0, k);
          VIEW (gsl_vector, view) b = FUNCTION (gsl_matrix, subrow) (dest, i, 0, k);

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
  else if (Uplo == CblasUpper)
    {
      for (i = 0; i < GSL_MIN(M, N - 1); i++)
        {
          VIEW (gsl_vector, const_view) a = FUNCTION (gsl_matrix, const_subrow) (src, i, i + 1, N - i - 1);
          VIEW (gsl_vector, view) b = FUNCTION (gsl_matrix, subrow) (dest, i, i + 1, N - i - 1);

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

    if (Uplo == CblasLower)
      {
        for (i = 1; i < M ; i++)
          {
            for (j = 0; j < GSL_MIN(i, N); j++)
              {
                for (k = 0; k < MULTIPLICITY; k++)
                  {
                    size_t e1 = (i *  dest_tda + j) * MULTIPLICITY + k ;
                    size_t e2 = (i *  src_tda + j) * MULTIPLICITY + k ;
                    dest->data[e1] = src->data[e2];
                  }
              }
          }
      }
    else if (Uplo == CblasUpper)
      {
        for (i = 0; i < M ; i++)
          {
            for (j = i + 1; j < N; j++)
              {
                for (k = 0; k < MULTIPLICITY; k++)
                  {
                    size_t e1 = (i *  dest_tda + j) * MULTIPLICITY + k ;
                    size_t e2 = (i *  src_tda + j) * MULTIPLICITY + k ;
                    dest->data[e1] = src->data[e2];
                  }
              }
          }
      }
    else
      {
        GSL_ERROR ("invalid Uplo parameter", GSL_EINVAL);
      }

    if (Diag == CblasNonUnit)
      {
        for (i = 0; i < GSL_MIN(M, N); i++)
          {
            for (k = 0; k < MULTIPLICITY; k++)
              {
                size_t e1 = (i *  dest_tda + i) * MULTIPLICITY + k ;
                size_t e2 = (i *  src_tda + i) * MULTIPLICITY + k ;
                dest->data[e1] = src->data[e2];
              }
          }
      }
  }

#endif

  return GSL_SUCCESS;
}


