/* spmatrix/prop_source.c
 * 
 * Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 Patrick Alken
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
FUNCTION (gsl_spmatrix, equal) (const TYPE (gsl_spmatrix) * a, const TYPE (gsl_spmatrix) * b)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (b->size1 != M || b->size2 != N)
    {
      GSL_ERROR_VAL("matrices must have same dimensions", GSL_EBADLEN, 0);
    }
  else if (a->sptype != b->sptype)
    {
      GSL_ERROR_VAL("trying to compare different sparse matrix types", GSL_EINVAL, 0);
    }
  else
    {
      const size_t nz = a->nz;
      size_t n, r;

      /* check for different number of non-zero elements */
      if (nz != b->nz)
        return 0;

      if (GSL_SPMATRIX_ISCOO(a))
        {
          /*
           * triplet formats could be out of order but identical, so use
           * gsl_spmatrix_get() on b for each aij
           */
          for (n = 0; n < nz; ++n)
            {
              ATOMIC * bptr = (ATOMIC *) FUNCTION (gsl_spmatrix, ptr) (b, a->i[n], a->p[n]);

              if (bptr == NULL)
                return 0;

              for (r = 0; r < MULTIPLICITY; ++r)
                {
                  if (a->data[MULTIPLICITY * n + r] != *(bptr + r))
                    return 0;
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSC(a))
        {
          /*
           * for CSC, both matrices should have everything
           * in the same order
           */

          /* check row indices and data */
          for (n = 0; n < nz; ++n)
            {
              if (a->i[n] != b->i[n])
                return 0;

              for (r = 0; r < MULTIPLICITY; ++r)
                {
                  if (a->data[MULTIPLICITY * n + r] != b->data[MULTIPLICITY * n + r])
                    return 0;
                }
            }

          /* check column pointers */
          for (n = 0; n < a->size2 + 1; ++n)
            {
              if (a->p[n] != b->p[n])
                return 0;
            }
        }
      else if (GSL_SPMATRIX_ISCSR(a))
        {
          /*
           * for CSR, both matrices should have everything
           * in the same order
           */

          /* check column indices and data */
          for (n = 0; n < nz; ++n)
            {
              if (a->i[n] != b->i[n])
                return 0;

              for (r = 0; r < MULTIPLICITY; ++r)
                {
                  if (a->data[MULTIPLICITY * n + r] != b->data[MULTIPLICITY * n + r])
                    return 0;
                }
            }

          /* check row pointers */
          for (n = 0; n < a->size1 + 1; ++n)
            {
              if (a->p[n] != b->p[n])
                return 0;
            }
        }
      else
        {
          GSL_ERROR_VAL("unknown sparse matrix type", GSL_EINVAL, 0);
        }

      return 1;
    }
}

#if !defined(UNSIGNED) && !defined(BASE_GSL_COMPLEX) && !defined(BASE_GSL_COMPLEX_FLOAT) && !defined(BASE_GSL_COMPLEX_LONG)

ATOMIC
FUNCTION (gsl_spmatrix, norm1) (const TYPE (gsl_spmatrix) * A)
{
  const size_t N = A->size2;
  ATOMIC value = (ATOMIC) 0;

  if (A->nz == 0)
    {
      return (ATOMIC) 0;
    }
  else if (GSL_SPMATRIX_ISCSC(A))
    {
      int * Ap = A->p;
      ATOMIC * Ad = A->data;
      size_t j;

      for (j = 0; j < N; ++j)
        {
          ATOMIC sum = (ATOMIC) 0;
          int p;

          for (p = Ap[j]; p < Ap[j + 1]; ++p)
            sum += (Ad[p] >= (ATOMIC) 0) ? Ad[p] : -Ad[p];

          if (sum > value)
            value = sum;
        }
    }
  else
    {
      ATOMIC * colsum = A->work.work_atomic;
      ATOMIC * Ad = A->data;
      size_t j;

      /* initialize column sums to zero */
      for (j = 0; j < N; ++j)
        colsum[j] = (ATOMIC) 0;

      if (GSL_SPMATRIX_ISCOO(A))
        {
          int * Ap = A->p;

          for (j = 0; j < A->nz; ++j)
            colsum[Ap[j]] += (Ad[j] >= (ATOMIC) 0) ? Ad[j] : -Ad[j];
        }
      else if (GSL_SPMATRIX_ISCSR(A))
        {
          int * Aj = A->i;

          for (j = 0; j < A->nz; ++j)
            colsum[Aj[j]] += (Ad[j] >= (ATOMIC) 0) ? Ad[j] : -Ad[j];
        }

      for (j = 0; j < N; ++j)
        {
          if (colsum[j] > value)
            value = colsum[j];
        }
    }

  return value;
}

#endif
