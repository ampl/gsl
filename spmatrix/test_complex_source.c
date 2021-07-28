/* spmatrix/test_complex_source.c
 * 
 * Copyright (C) 2018, 2019, 2020 Patrick Alken
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
test_random()
  Create a random sparse matrix with approximately
M*N*density non-zero entries in [lower,upper]

Inputs: M       - number of rows
        N       - number of columns
        density - sparse density \in [0,1]
                  0 = no non-zero entries
                  1 = all m*n entries are filled
        lower   - lower bound on entries
        upper   - upper bound on entries
        r       - random number generator

Return: pointer to sparse matrix in triplet format (must be freed by caller)

Notes:
1) non-zero matrix entries are uniformly distributed in [0,1]
*/

static TYPE (gsl_spmatrix) *
FUNCTION (test, random)(const size_t M, const size_t N, const double density,
                        const double lower, const double upper, const gsl_rng *r)
{
  size_t nnzwanted = (size_t) floor(M * N * GSL_MIN(density, 1.0));
  TYPE (gsl_spmatrix) * m = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, nnzwanted, GSL_SPMATRIX_COO);

  while (FUNCTION (gsl_spmatrix, nnz) (m) < nnzwanted)
    {
      /* generate a random row and column */
      size_t i = gsl_rng_uniform(r) * M;
      size_t j = gsl_rng_uniform(r) * N;

      /* generate random m_{ij} and add it */
      BASE z;
      GSL_REAL (z) = (ATOMIC) (gsl_rng_uniform(r) * (upper - lower) + lower);
      GSL_IMAG (z) = (ATOMIC) (gsl_rng_uniform(r) * (upper - lower) + lower);
      FUNCTION (gsl_spmatrix, set) (m, i, j, z);
    }

  return m;
}

static TYPE (gsl_spmatrix) *
FUNCTION (test, random_int)(const size_t M, const size_t N, const double density,
                            const double lower, const double upper, const gsl_rng *r)
{
  size_t nnzwanted = (size_t) floor(M * N * GSL_MIN(density, 1.0));
  TYPE (gsl_spmatrix) * m = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, nnzwanted, GSL_SPMATRIX_COO);

  while (FUNCTION (gsl_spmatrix, nnz) (m) < nnzwanted)
    {
      /* generate a random row and column */
      size_t i = gsl_rng_uniform(r) * M;
      size_t j = gsl_rng_uniform(r) * N;

      /* generate random m_{ij} and add it */
      int xr = (int) (gsl_rng_uniform(r) * (upper - lower) + lower);
      int xi = (int) (gsl_rng_uniform(r) * (upper - lower) + lower);
      BASE x;

      GSL_REAL(x) = (ATOMIC) xr;
      GSL_IMAG(x) = (ATOMIC) xi;
      FUNCTION (gsl_spmatrix, set) (m, i, j, x);
    }

  return m;
}

static void
FUNCTION (test, random_dense)(TYPE (gsl_matrix) * m, const double lower,
                              const double upper, const gsl_rng *r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double x;
          BASE z;
          
          x = gsl_rng_uniform(r) * (upper - lower) + lower;
          GSL_REAL (z) = (ATOMIC) x;

          x = gsl_rng_uniform(r) * (upper - lower) + lower;
          GSL_IMAG (z) = (ATOMIC) x;

          FUNCTION (gsl_matrix, set) (m, i, j, z);
        }
    }
}

static void
FUNCTION (test, random_vector)(TYPE (gsl_vector) * x, const double lower,
                               const double upper, const gsl_rng *r)
{
  size_t i;

  for (i = 0; i < x->size; ++i)
    {
      BASE z;
      GSL_REAL (z) = gsl_rng_uniform(r) * (upper - lower) + lower;
      GSL_IMAG (z) = gsl_rng_uniform(r) * (upper - lower) + lower;
      FUNCTION (gsl_vector, set) (x, i, z);
    }
}

static void
FUNCTION (test, alloc) (const size_t M, const size_t N, const int sptype)
{
  TYPE (gsl_spmatrix) * m = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, 10, sptype);
  size_t i, j;

  gsl_test (m->data == 0, NAME (gsl_spmatrix) "_alloc(%s) returns valid data pointer", FUNCTION (gsl_spmatrix, type) (m));
  gsl_test (m->i == 0, NAME (gsl_spmatrix) "_alloc(%s) returns valid i pointer", FUNCTION (gsl_spmatrix, type) (m));
  gsl_test (m->p == 0, NAME (gsl_spmatrix) "_alloc(%s) returns valid p pointer", FUNCTION (gsl_spmatrix, type) (m));
  gsl_test (m->size1 != M, NAME (gsl_spmatrix) "_alloc(%s) returns valid size1", FUNCTION (gsl_spmatrix, type) (m));
  gsl_test (m->size2 != N, NAME (gsl_spmatrix) "_alloc(%s) returns valid size2", FUNCTION (gsl_spmatrix, type) (m));
  gsl_test (m->nz != 0, NAME (gsl_spmatrix) "_alloc(%s) returns valid nz", FUNCTION (gsl_spmatrix, type) (m));
  gsl_test (m->sptype != sptype, NAME (gsl_spmatrix) "_alloc(%s) returns valid sptype", FUNCTION (gsl_spmatrix, type) (m));

  /* test that a newly allocated matrix is initialized to zero */
  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE mij = FUNCTION (gsl_spmatrix, get) (m, i, j);

          if (GSL_REAL(mij) != (ATOMIC) 0 || GSL_IMAG(mij) != (ATOMIC) 0)
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_alloc(%s) initializes matrix to zero", FUNCTION (gsl_spmatrix, type) (m));

  FUNCTION (gsl_spmatrix, free) (m);
}

static void
FUNCTION (test, realloc) (const size_t M, const size_t N)
{
  TYPE (gsl_spmatrix) * m = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, 10, GSL_SPMATRIX_COO);
  const size_t nnz_total = 30;
  BASE z;
  size_t i, j, nnz;

  GSL_REAL (z) = (ATOMIC) 5;
  GSL_IMAG (z) = (ATOMIC) 200;

  /* add enough elements to m to trigger a realloc call, then check the result */

  nnz = 0;
  i = 0;
  j = 0;
  while (nnz++ < nnz_total)
    {
      if (i >= m->size1)
        {
          i = 0;
          if (++j >= m->size2)
            break;
        }

      FUNCTION (gsl_spmatrix, set) (m, i++, j, z);
    }

  status = 0;
  nnz = 0;
  i = 0;
  j = 0;
  while (nnz++ < nnz_total)
    {
      BASE y;

      if (i >= m->size1)
        {
          i = 0;
          if (++j >= m->size2)
            break;
        }

      y = FUNCTION (gsl_spmatrix, get) (m, i++, j);

      if (GSL_REAL(y) != GSL_REAL(z) || GSL_IMAG(y) != GSL_IMAG(z))
        status = 1;
    }

  gsl_test (status, NAME (gsl_spmatrix) "_realloc[%zu,%zu]", M, N);

  FUNCTION (gsl_spmatrix, free) (m);
}

static void
FUNCTION (test, getset) (const size_t M, const size_t N, const int sptype)
{
  TYPE (gsl_spmatrix) * A = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, M * N, GSL_SPMATRIX_COO);
  TYPE (gsl_spmatrix) * B = NULL;
  size_t i, j;
  size_t k = 0;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE z;
          GSL_REAL(z) = (ATOMIC) k++;
          GSL_IMAG(z) = (ATOMIC) (k + 1000);
          FUNCTION (gsl_spmatrix, set) (A, i, j, z);
        }
    }

  status = 0;
  k = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE z;
          GSL_REAL(z) = (ATOMIC) k++;
          GSL_IMAG(z) = (ATOMIC) (k + 1000);
          if ((A->i[k - 1] != (int) i) || (A->p[k - 1] != (int) j) ||
              (A->data[2*(k - 1)] != GSL_REAL(z)) || (A->data[2*(k - 1) + 1] != GSL_IMAG(z)))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_set[%zu,%zu] writes into arrays", M, N);

  B = FUNCTION (gsl_spmatrix, compress) (A, sptype);

  status = 0;
  k = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE y = FUNCTION (gsl_spmatrix, get) (B, i, j);
          BASE z;
          GSL_REAL(z) = (ATOMIC) k++;
          GSL_IMAG(z) = (ATOMIC) (k + 1000);

          if (GSL_REAL(y) != GSL_REAL(z) || GSL_IMAG(y) != GSL_IMAG(z))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_get[%zu,%zu](%s) reads from arrays",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  status = FUNCTION (gsl_spmatrix, nnz) (B) != M * N;
  gsl_test (status, NAME (gsl_spmatrix) "_nnz[%zu,%zu](%s) %zu",
            M, N, FUNCTION (gsl_spmatrix, type) (B), FUNCTION (gsl_spmatrix, nnz) (B));

  /* test gsl_spmatrix_ptr */
  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE mij = FUNCTION (gsl_spmatrix, get) (B, i, j);
          BASE * ptr = FUNCTION (gsl_spmatrix, ptr) (B, i, j);
          BASE mij2;

          GSL_REAL(*ptr) += (ATOMIC) 2;
          GSL_IMAG(*ptr) += (ATOMIC) 5;
          mij2 = FUNCTION (gsl_spmatrix, get) (B, i, j);

          if ((GSL_REAL(mij2) != GSL_REAL(mij) + (ATOMIC) 2) ||
              (GSL_IMAG(mij2) != GSL_IMAG(mij) + (ATOMIC) 5))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_ptr[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  /* test set_zero */
  FUNCTION (gsl_spmatrix, set_zero) (B);
  status = FUNCTION (gsl_spmatrix, nnz) (B) != 0;
  gsl_test (status, NAME (gsl_spmatrix) "_set_zero[%zu,%zu](%s) nnz=%zu",
            M, N, FUNCTION (gsl_spmatrix, type) (B), FUNCTION (gsl_spmatrix, nnz) (B));

  if (B->sptype == GSL_SPMATRIX_COO)
    {
      const size_t min = GSL_MIN(M, N);
      const size_t expected_nnz = min;
      BASE one = ONE;
      BASE zero = ZERO;
      size_t nnz;
      BASE z;

      /* test setting duplicate values */
      status = 0;
      k = 0;
      for (i = 0; i < min; ++i)
        {
          for (j = 0; j < 5; ++j)
            {
              BASE x, y;

              GSL_REAL(x) = (ATOMIC) k++;
              GSL_IMAG(x) = (ATOMIC) (k + 1000);

              FUNCTION (gsl_spmatrix, set) (B, i, i, x);
              y = FUNCTION (gsl_spmatrix, get) (B, i, i);

              if (GSL_REAL(x) != GSL_REAL(y) || GSL_IMAG(x) != GSL_IMAG(y))
                status = 1;
            }
        }

      gsl_test (status, NAME (gsl_spmatrix) " duplicates[%zu,%zu]", M, N);

      nnz = FUNCTION (gsl_spmatrix, nnz) (B);
      status = nnz != expected_nnz;
      gsl_test (status, NAME (gsl_spmatrix) " duplicates nnz[%zu,%zu] %zu/%zu",
                M, N, nnz, expected_nnz);

      /* test setting an element to 0 */
      FUNCTION (gsl_spmatrix, set) (B, 0, 0, one);
      FUNCTION (gsl_spmatrix, set) (B, 0, 0, zero);
      z = FUNCTION (gsl_spmatrix, get) (B, 0, 0);
      status = GSL_REAL(z) != (ATOMIC) 0 || GSL_IMAG(z) != (ATOMIC) 0;
      gsl_test (status, NAME (gsl_spmatrix) " zero element[%zu,%zu] %f",
                M, N, FUNCTION (gsl_spmatrix, get) (B, 0, 0));
    }

  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
}

static void
FUNCTION (test, memcpy) (const size_t M, const size_t N, const int sptype,
                         const double density, gsl_rng * r)
{
  /* test memcpy */
  {
    TYPE (gsl_spmatrix) * m = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
    TYPE (gsl_spmatrix) * A = FUNCTION (gsl_spmatrix, compress) (m, sptype);
    TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, m->nz, sptype);

    FUNCTION (gsl_spmatrix, memcpy) (B, A);

    status = FUNCTION (gsl_spmatrix, equal) (A, B) != 1;
    gsl_test (status, NAME (gsl_spmatrix) "_memcpy[%zu,%zu](%s) equality",
              M, N, FUNCTION (gsl_spmatrix, type) (A));

    FUNCTION (gsl_spmatrix, free) (m);
    FUNCTION (gsl_spmatrix, free) (A);
    FUNCTION (gsl_spmatrix, free) (B);
  }

  /* test transpose memcpy */
  {
    TYPE (gsl_spmatrix) * m = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
    TYPE (gsl_spmatrix) * A = FUNCTION (gsl_spmatrix, compress) (m, sptype);
    TYPE (gsl_spmatrix) * AT = FUNCTION (gsl_spmatrix, alloc_nzmax) (N, M, m->nz, sptype);
    size_t i, j;

    FUNCTION (gsl_spmatrix, transpose_memcpy) (AT, A);

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            BASE Aij = FUNCTION (gsl_spmatrix, get) (A, i, j);
            BASE ATji = FUNCTION (gsl_spmatrix, get) (AT, j, i);

            if (GSL_REAL(Aij) != GSL_REAL(ATji) ||
                GSL_IMAG(Aij) != GSL_IMAG(ATji))
              status = 1;
          }
      }

    gsl_test(status, NAME (gsl_spmatrix) "_transpose_memcpy[%zu,%zu](%s)",
             M, N, FUNCTION (gsl_spmatrix, type) (A));

    FUNCTION (gsl_spmatrix, free) (m);
    FUNCTION (gsl_spmatrix, free) (A);
    FUNCTION (gsl_spmatrix, free) (AT);
  }
}

static void
FUNCTION (test, transpose) (const size_t M, const size_t N, const int sptype,
                            const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * A = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (A, sptype);
  size_t i;

  FUNCTION (gsl_spmatrix, transpose) (B);

  status = 0;
  for (i = 0; i < A->nz; ++i)
    {
      BASE Aij = *(BASE *) &A->data[2 * i];
      BASE Bji = FUNCTION (gsl_spmatrix, get) (B, A->p[i], A->i[i]);

      if (GSL_REAL(Aij) != GSL_REAL(Bji) ||
          GSL_IMAG(Aij) != GSL_IMAG(Bji))
        status = 1;
    }

  gsl_test(status, NAME (gsl_spmatrix) "_transpose[%zu,%zu](%s)",
           M, N, FUNCTION (gsl_spmatrix, type) (B));

  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
}

static void
FUNCTION (test, scale) (const size_t M, const size_t N, const int sptype,
                        const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * A = FUNCTION (test, random_int) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (A, sptype);
  TYPE (gsl_vector) * x;
  BASE s;
  size_t i, j;

  GSL_SET_COMPLEX (&s, 2.0, 3.0);
  FUNCTION (gsl_spmatrix, scale) (B, s);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE aij = FUNCTION (gsl_spmatrix, get) (A, i, j);
          BASE bij = FUNCTION (gsl_spmatrix, get) (B, i, j);
          BASE cij;

          GSL_SET_COMPLEX(&cij,
                          GSL_REAL(aij)*GSL_REAL(s) - GSL_IMAG(aij)*GSL_IMAG(s),
                          GSL_IMAG(aij)*GSL_REAL(s) + GSL_REAL(aij)*GSL_IMAG(s));

          if (GSL_REAL(bij) != GSL_REAL(cij) ||
              GSL_IMAG(bij) != GSL_IMAG(cij))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_scale[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  /* reset B = A */
  FUNCTION (gsl_spmatrix, free) (B);
  B = FUNCTION (gsl_spmatrix, compress) (A, sptype);

  /* test column scaling */
  x = FUNCTION (gsl_vector, alloc) (A->size2);
  FUNCTION (test, random_vector) (x, 1.0, 20.0, r);

  FUNCTION (gsl_spmatrix, scale_columns) (B, x);

  status = 0;
  for (i = 0; i < A->nz; ++i)
    {
      BASE aij = *(BASE *) &A->data[2 * i];
      BASE bij = FUNCTION (gsl_spmatrix, get) (B, A->i[i], A->p[i]);
      BASE xj = FUNCTION (gsl_vector, get) (x, A->p[i]);
      BASE cij; /* aij * xj */

      GSL_SET_COMPLEX(&cij,
                      GSL_REAL(aij)*GSL_REAL(xj) - GSL_IMAG(aij)*GSL_IMAG(xj),
                      GSL_IMAG(aij)*GSL_REAL(xj) + GSL_REAL(aij)*GSL_IMAG(xj));

      gsl_test_rel(GSL_REAL(bij), GSL_REAL(cij), 10.0 * GSL_DBL_EPSILON,
                   NAME (gsl_spmatrix) "_scale_columns[%zu,%zu](%s) real",
                   M, N, FUNCTION (gsl_spmatrix, type) (B));

      gsl_test_rel(GSL_IMAG(bij), GSL_IMAG(cij), 10.0 * GSL_DBL_EPSILON,
                   NAME (gsl_spmatrix) "_scale_columns[%zu,%zu](%s) imag",
                   M, N, FUNCTION (gsl_spmatrix, type) (B));
    }

  FUNCTION (gsl_vector, free) (x);

  /* reset B = A */
  FUNCTION (gsl_spmatrix, free) (B);
  B = FUNCTION (gsl_spmatrix, compress) (A, sptype);

  /* test row scaling */
  x = FUNCTION (gsl_vector, alloc) (A->size1);
  FUNCTION (test, random_vector) (x, 1.0, 20.0, r);

  FUNCTION (gsl_spmatrix, scale_rows) (B, x);

  status = 0;
  for (i = 0; i < A->nz; ++i)
    {
      BASE aij = *(BASE *) &A->data[2 * i];
      BASE bij = FUNCTION (gsl_spmatrix, get) (B, A->i[i], A->p[i]);
      BASE xi = FUNCTION (gsl_vector, get) (x, A->i[i]);
      BASE cij; /* aij * xi */

      GSL_SET_COMPLEX(&cij,
                      GSL_REAL(aij)*GSL_REAL(xi) - GSL_IMAG(aij)*GSL_IMAG(xi),
                      GSL_IMAG(aij)*GSL_REAL(xi) + GSL_REAL(aij)*GSL_IMAG(xi));

      gsl_test_rel(GSL_REAL(bij), GSL_REAL(cij), 10.0 * GSL_DBL_EPSILON,
                   NAME (gsl_spmatrix) "_scale_rows[%zu,%zu](%s) real",
                   M, N, FUNCTION (gsl_spmatrix, type) (B));

      gsl_test_rel(GSL_IMAG(bij), GSL_IMAG(cij), 10.0 * GSL_DBL_EPSILON,
                   NAME (gsl_spmatrix) "_scale_rows[%zu,%zu](%s) imag",
                   M, N, FUNCTION (gsl_spmatrix, type) (B));
    }

  FUNCTION (gsl_vector, free) (x);
  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
}

static void
FUNCTION (test, add) (const size_t M, const size_t N, const int sptype,
                      const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * m1 = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * m2 = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * A = FUNCTION (gsl_spmatrix, compress) (m1, sptype);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (m2, sptype);
  TYPE (gsl_spmatrix) * C = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, 1, sptype);
  size_t i, j;

  FUNCTION (gsl_spmatrix, add) (C, A, B);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE aij = FUNCTION (gsl_spmatrix, get) (A, i, j);
          BASE bij = FUNCTION (gsl_spmatrix, get) (B, i, j);
          BASE cij = FUNCTION (gsl_spmatrix, get) (C, i, j);

          if ((GSL_REAL(aij) + GSL_REAL(bij) != GSL_REAL(cij)) ||
              (GSL_IMAG(aij) + GSL_IMAG(bij) != GSL_IMAG(cij)))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_add[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (A));

  /* test again with C = 2*A */
  FUNCTION (gsl_spmatrix, add) (C, A, A);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE aij = FUNCTION (gsl_spmatrix, get) (A, i, j);
          BASE cij = FUNCTION (gsl_spmatrix, get) (C, i, j);

          if ((2 * GSL_REAL(aij) != GSL_REAL(cij)) ||
              (2 * GSL_IMAG(aij) != GSL_IMAG(cij)))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_add[%zu,%zu](%s) duplicate",
            M, N, FUNCTION (gsl_spmatrix, type) (A));

  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
  FUNCTION (gsl_spmatrix, free) (C);
  FUNCTION (gsl_spmatrix, free) (m1);
  FUNCTION (gsl_spmatrix, free) (m2);
}

static void
FUNCTION (test, dense_add) (const size_t M, const size_t N, const int sptype,
                            const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * m = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (m, sptype);
  TYPE (gsl_matrix) * A = FUNCTION (gsl_matrix, alloc) (M, N);
  TYPE (gsl_matrix) * A_copy = FUNCTION (gsl_matrix, alloc) (M, N);
  size_t i, j;

  FUNCTION (test, random_dense) (A, 1.0, 20.0, r);
  FUNCTION (gsl_matrix, memcpy) (A_copy, A);

  FUNCTION (gsl_spmatrix, dense_add) (A, B);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE aij = FUNCTION (gsl_matrix, get) (A_copy, i, j);
          BASE bij = FUNCTION (gsl_spmatrix, get) (B, i, j);
          BASE cij = FUNCTION (gsl_matrix, get) (A, i, j);

          if ((GSL_REAL(aij) + GSL_REAL(bij) != GSL_REAL(cij)) ||
              (GSL_IMAG(aij) + GSL_IMAG(bij) != GSL_IMAG(cij)))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_dense_add[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  FUNCTION (gsl_matrix, free) (A);
  FUNCTION (gsl_matrix, free) (A_copy);
  FUNCTION (gsl_spmatrix, free) (B);
  FUNCTION (gsl_spmatrix, free) (m);
}

static void
FUNCTION (test, dense_sub) (const size_t M, const size_t N, const int sptype,
                            const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * m = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (m, sptype);
  TYPE (gsl_matrix) * A = FUNCTION (gsl_matrix, alloc) (M, N);
  TYPE (gsl_matrix) * A_copy = FUNCTION (gsl_matrix, alloc) (M, N);
  size_t i, j;

  FUNCTION (test, random_dense) (A, 1.0, 20.0, r);
  FUNCTION (gsl_matrix, memcpy) (A_copy, A);

  FUNCTION (gsl_spmatrix, dense_sub) (A, B);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE aij = FUNCTION (gsl_matrix, get) (A_copy, i, j);
          BASE bij = FUNCTION (gsl_spmatrix, get) (B, i, j);
          BASE cij = FUNCTION (gsl_matrix, get) (A, i, j);

          if ((GSL_REAL(aij) - GSL_REAL(bij) != GSL_REAL(cij)) ||
              (GSL_IMAG(aij) - GSL_IMAG(bij) != GSL_IMAG(cij)))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_dense_sub[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  FUNCTION (gsl_matrix, free) (A);
  FUNCTION (gsl_matrix, free) (A_copy);
  FUNCTION (gsl_spmatrix, free) (B);
  FUNCTION (gsl_spmatrix, free) (m);
}

static void
FUNCTION (test, convert) (const size_t M, const size_t N, const int sptype,
                          const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * A = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (A, sptype);
  TYPE (gsl_spmatrix) * C = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, 1, GSL_SPMATRIX_COO);
  TYPE (gsl_matrix) * D = FUNCTION (gsl_matrix, alloc) (M, N);
  size_t i, j;

  /* convert D := B */
  FUNCTION (gsl_spmatrix, sp2d) (D, B);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE bij = FUNCTION (gsl_spmatrix, get) (B, i, j);
          BASE dij = FUNCTION (gsl_matrix, get) (D, i, j);

          if ((GSL_REAL(bij) != GSL_REAL(dij)) ||
              (GSL_IMAG(bij) != GSL_IMAG(dij)))
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_sp2d[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  /* convert C := D */
  FUNCTION (gsl_spmatrix, d2sp) (C, D);

  status = FUNCTION (gsl_spmatrix, equal) (A, C) != 1;
  gsl_test (status, NAME (gsl_spmatrix) "_d2sp[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
  FUNCTION (gsl_spmatrix, free) (C);
  FUNCTION (gsl_matrix, free) (D);
}

static void
FUNCTION (test, io_ascii) (const size_t M, const size_t N, const int sptype,
                           const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * A = FUNCTION (test, random_int) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (A, sptype);
  TYPE (gsl_spmatrix) * C;
  char filename[] = "test.dat";
  FILE *f;

  f = fopen (filename, "w");
  FUNCTION (gsl_spmatrix, fprintf) (f, B, OUT_FORMAT);
  fclose (f);

  f = fopen (filename, "r");
  C = FUNCTION (gsl_spmatrix, fscanf) (f);
  fclose (f);

  status = FUNCTION (gsl_spmatrix, equal) (A, C) != 1;
  gsl_test (status, NAME (gsl_spmatrix) "_fscanf[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  unlink (filename);

  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
  FUNCTION (gsl_spmatrix, free) (C);
}

static void
FUNCTION (test, io_binary) (const size_t M, const size_t N, const int sptype,
                            const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * A = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (A, sptype);
  TYPE (gsl_spmatrix) * C = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, A->nz, sptype);
  char filename[] = "test.dat";
  FILE *f;

  f = fopen (filename, "wb");
  FUNCTION (gsl_spmatrix, fwrite) (f, B);
  fclose (f);

  f = fopen (filename, "rb");
  FUNCTION (gsl_spmatrix, fread) (f, C);
  fclose (f);

  status = FUNCTION (gsl_spmatrix, equal) (B, C) != 1;
  gsl_test (status, NAME (gsl_spmatrix) "_fread[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  unlink (filename);

  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
  FUNCTION (gsl_spmatrix, free) (C);
}

static void
FUNCTION (test, all) (const size_t M, const size_t N, const double density, gsl_rng * r)
{
  FUNCTION (test, alloc) (M, N, GSL_SPMATRIX_COO);
  FUNCTION (test, alloc) (M, N, GSL_SPMATRIX_CSC);
  FUNCTION (test, alloc) (M, N, GSL_SPMATRIX_CSR);

  FUNCTION (test, realloc) (M, N);

  FUNCTION (test, getset) (M, N, GSL_SPMATRIX_COO);
  FUNCTION (test, getset) (M, N, GSL_SPMATRIX_CSC);
  FUNCTION (test, getset) (M, N, GSL_SPMATRIX_CSR);

  FUNCTION (test, memcpy) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, memcpy) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, memcpy) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, transpose) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, transpose) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, transpose) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, scale) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, scale) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, scale) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, add) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, add) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, dense_add) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, dense_add) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, dense_add) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, dense_sub) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, dense_sub) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, dense_sub) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, convert) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, convert) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, convert) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, io_ascii) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, io_ascii) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, io_ascii) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, io_binary) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, io_binary) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, io_binary) (M, N, GSL_SPMATRIX_CSR, density, r);
}
