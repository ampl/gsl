/* matrix/prop_source.c
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
FUNCTION (gsl_matrix, equal) (const TYPE (gsl_matrix) * a, const TYPE (gsl_matrix) * b)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (b->size1 != M || b->size2 != N)
    {
      GSL_ERROR_VAL ("matrices must have same dimensions", GSL_EBADLEN, 0);
    }
  else 
    {
      const size_t tda_a = a->tda;
      const size_t tda_b = b->tda;

      size_t i, j, k;

      for (i = 0; i < M; i++)
        {
          for (j = 0; j < N; j++)
            {
              for (k = 0; k < MULTIPLICITY; k++) 
                {
                  if (a->data[(i * tda_a + j) * MULTIPLICITY + k] 
                      != b->data[(i * tda_b + j) * MULTIPLICITY + k])
                    {
                      return 0;
                    }
                }
            }
        }
    }
  return 1;
}


int
FUNCTION (gsl_matrix, isnull) (const TYPE (gsl_matrix) * m)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;
  const size_t tda = m->tda ;
  
  size_t i, j, k;

  for (i = 0; i < size1 ; i++)
    {
      for (j = 0; j < size2; j++)
        {
          for (k = 0; k < MULTIPLICITY; k++) 
            {
              if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
                {
                  return 0;
                }
            }
        }
    }
      
  return 1;
}


int
FUNCTION (gsl_matrix, ispos) (const TYPE (gsl_matrix) * m)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;
  const size_t tda = m->tda ;
  
  size_t i, j, k;

  for (i = 0; i < size1 ; i++)
    {
      for (j = 0; j < size2; j++)
        {
          for (k = 0; k < MULTIPLICITY; k++) 
            {
              if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
                {
                  return 0;
                }
            }
        }
    }
      
  return 1;
}


int
FUNCTION (gsl_matrix, isneg) (const TYPE (gsl_matrix) * m)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;
  const size_t tda = m->tda ;
  
  size_t i, j, k;

  for (i = 0; i < size1 ; i++)
    {
      for (j = 0; j < size2; j++)
        {
          for (k = 0; k < MULTIPLICITY; k++) 
            {
              if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
                {
                  return 0;
                }
            }
        }
    }
      
  return 1;
}


int
FUNCTION (gsl_matrix, isnonneg) (const TYPE (gsl_matrix) * m)
{
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;
  const size_t tda = m->tda ;
  
  size_t i, j, k;

  for (i = 0; i < size1 ; i++)
    {
      for (j = 0; j < size2; j++)
        {
          for (k = 0; k < MULTIPLICITY; k++) 
            {
              if (m->data[(i * tda + j) * MULTIPLICITY + k] < 0.0)
                {
                  return 0;
                }
            }
        }
    }
      
  return 1;
}

#if !defined(UNSIGNED) && !defined(BASE_GSL_COMPLEX) && !defined(BASE_GSL_COMPLEX_FLOAT) && !defined(BASE_GSL_COMPLEX_LONG)

ATOMIC
FUNCTION (gsl_matrix, norm1) (const TYPE (gsl_matrix) * m)
{
  ATOMIC value = (ATOMIC) 0;
  size_t j;

  for (j = 0; j < m->size2; ++j)
    {
      VIEW (gsl_vector, const_view) mj = FUNCTION (gsl_matrix, const_column) (m, j);
      ATOMIC sum;
      
#if defined(BASE_DOUBLE)
      sum = gsl_blas_dasum(&mj.vector);
#elif defined(BASE_FLOAT)
      sum = gsl_blas_sasum(&mj.vector);
#else
      {
        size_t i;

        sum = (ATOMIC) 0;

        for (i = 0; i < m->size1; ++i)
          {
            ATOMIC mij = FUNCTION (gsl_vector, get) (&mj.vector, i);

            if (mij >= (ATOMIC) 0)
              sum += mij;
            else
              sum += -mij;
          }
      }
#endif

      if (sum > value)
        value = sum;
    }

  return value;
}

#endif
