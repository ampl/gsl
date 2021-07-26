/* spmatrix/test_source.c
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
      BASE x = (BASE) (gsl_rng_uniform(r) * (upper - lower) + lower);
      FUNCTION (gsl_spmatrix, set) (m, i, j, x);
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
      int x = (int) (gsl_rng_uniform(r) * (upper - lower) + lower);
      FUNCTION (gsl_spmatrix, set) (m, i, j, (BASE) x);
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
          double mij = gsl_rng_uniform(r) * (upper - lower) + lower;
          FUNCTION (gsl_matrix, set) (m, i, j, (BASE) mij);
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
      double xi = gsl_rng_uniform(r) * (upper - lower) + lower;
      FUNCTION (gsl_vector, set) (x, i, (BASE) xi);
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

          if (mij != (BASE) 0)
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
  const BASE x = (BASE) 5;
  const size_t nnz_total = 30;
  size_t i, j, nnz;

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

      FUNCTION (gsl_spmatrix, set) (m, i++, j, x);
    }

  status = 0;
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

      if (FUNCTION (gsl_spmatrix, get) (m, i++, j) != x)
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
  size_t i, j, idx;
  size_t k = 0;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          k++;
          FUNCTION (gsl_spmatrix, set) (A, i, j, (BASE) k);
        }
    }

  status = 0;
  k = 0;
  idx = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          k++;
          if ((A->i[idx] != (int) i) || (A->p[idx] != (int) j) ||
              (A->data[idx] != (BASE) k))
            status = 1;
          idx++;
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
          k++;
          if (FUNCTION (gsl_spmatrix, get) (B, i, j) != (BASE) k)
            status = 1;
        }
    }

  gsl_test (status, NAME (gsl_spmatrix) "_get[%zu,%zu](%s) reads from arrays",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  status = FUNCTION (gsl_spmatrix, nnz) (B) != M * N;
  gsl_test (status, NAME (gsl_spmatrix) "_nnz[%zu,%zu](%s) %zu/%zu",
            M, N, FUNCTION (gsl_spmatrix, type) (B),
            FUNCTION (gsl_spmatrix, nnz) (B), M * N);

  /* test gsl_spmatrix_ptr */
  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE mij = FUNCTION (gsl_spmatrix, get) (B, i, j);
          BASE * ptr = FUNCTION (gsl_spmatrix, ptr) (B, i, j);

          if (ptr != NULL)
            {
              *ptr += (ATOMIC) 2;
              if (FUNCTION (gsl_spmatrix, get) (B, i, j) != (BASE) (mij + (ATOMIC) 2))
                status = 1;
            }
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
      size_t nnz;

      /* test setting duplicate values */
      status = 0;
      k = 0;
      for (i = 0; i < min; ++i)
        {
          for (j = 0; j < 5; ++j)
            {
              BASE x = (BASE) ++k;
              BASE y;

              FUNCTION (gsl_spmatrix, set) (B, i, i, x);
              y = FUNCTION (gsl_spmatrix, get) (B, i, i);
              if (x != y)
                status = 1;
            }
        }

      gsl_test (status, NAME (gsl_spmatrix) " duplicates[%zu,%zu]", M, N);

      nnz = FUNCTION (gsl_spmatrix, nnz) (B);
      status = nnz != expected_nnz;
      gsl_test (status, NAME (gsl_spmatrix) " duplicates nnz[%zu,%zu] %zu/%zu",
                M, N, nnz, expected_nnz);

      /* test setting an element to 0 */
      FUNCTION (gsl_spmatrix, set) (B, 0, 0, (BASE) 1);
      FUNCTION (gsl_spmatrix, set) (B, 0, 0, (BASE) 0);
      status = FUNCTION (gsl_spmatrix, get) (B, 0, 0) != (BASE) 0;
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
    TYPE (gsl_spmatrix) * C = FUNCTION (gsl_spmatrix, alloc_nzmax) (M, N, 1, sptype);

    FUNCTION (gsl_spmatrix, memcpy) (B, A);
    status = FUNCTION (gsl_spmatrix, equal) (A, B) != 1;
    gsl_test (status, NAME (gsl_spmatrix) "_memcpy[%zu,%zu](%s) equality",
              M, N, FUNCTION (gsl_spmatrix, type) (A));

    /* this will force a realloc call */
    FUNCTION (gsl_spmatrix, memcpy) (C, A);
    status = FUNCTION (gsl_spmatrix, equal) (A, C) != 1;
    gsl_test (status, NAME (gsl_spmatrix) "_memcpy[%zu,%zu](%s) realloc equality",
              M, N, FUNCTION (gsl_spmatrix, type) (A));

    FUNCTION (gsl_spmatrix, free) (m);
    FUNCTION (gsl_spmatrix, free) (A);
    FUNCTION (gsl_spmatrix, free) (B);
    FUNCTION (gsl_spmatrix, free) (C);
  }

  /* test transpose memcpy */
  {
    TYPE (gsl_spmatrix) * m = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
    TYPE (gsl_spmatrix) * A = FUNCTION (gsl_spmatrix, compress) (m, sptype);
    TYPE (gsl_spmatrix) * AT = FUNCTION (gsl_spmatrix, alloc_nzmax) (N, M, m->nz, sptype);
    TYPE (gsl_spmatrix) * BT = FUNCTION (gsl_spmatrix, alloc_nzmax) (N, M, 1, sptype);
    size_t i, j;

    FUNCTION (gsl_spmatrix, transpose_memcpy) (AT, A);
    FUNCTION (gsl_spmatrix, transpose_memcpy) (BT, A); /* force realloc call */

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            BASE Aij = FUNCTION (gsl_spmatrix, get) (A, i, j);
            BASE ATji = FUNCTION (gsl_spmatrix, get) (AT, j, i);
            BASE BTji = FUNCTION (gsl_spmatrix, get) (BT, j, i);

            if (Aij != ATji)
              status = 1;
            if (Aij != BTji)
              status = 2;
          }
      }

    gsl_test(status == 1, NAME (gsl_spmatrix) "_transpose_memcpy[%zu,%zu](%s) AT",
             M, N, FUNCTION (gsl_spmatrix, type) (A));
    gsl_test(status == 2, NAME (gsl_spmatrix) "_transpose_memcpy[%zu,%zu](%s) BT",
             M, N, FUNCTION (gsl_spmatrix, type) (A));

    FUNCTION (gsl_spmatrix, free) (m);
    FUNCTION (gsl_spmatrix, free) (A);
    FUNCTION (gsl_spmatrix, free) (AT);
    FUNCTION (gsl_spmatrix, free) (BT);
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
      BASE Aij = A->data[i];
      BASE Bji = FUNCTION (gsl_spmatrix, get) (B, A->p[i], A->i[i]);

      if (Aij != Bji)
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
  TYPE (gsl_spmatrix) * A = FUNCTION (test, random) (M, N, density, 1.0, 20.0, r);
  TYPE (gsl_spmatrix) * B = FUNCTION (gsl_spmatrix, compress) (A, sptype);
  TYPE (gsl_vector) * x;
  size_t i, j;

  FUNCTION (gsl_spmatrix, scale) (B, 2.0);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          BASE aij = FUNCTION (gsl_spmatrix, get) (A, i, j);
          BASE bij = FUNCTION (gsl_spmatrix, get) (B, i, j);
          if (bij != (ATOMIC) (2*aij))
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
      BASE aij = A->data[i];
      BASE bij = FUNCTION (gsl_spmatrix, get) (B, A->i[i], A->p[i]);
      BASE xj = FUNCTION (gsl_vector, get) (x, A->p[i]);

      if (bij != (ATOMIC) (xj * aij))
        status = 1;
    }

  gsl_test (status, NAME (gsl_spmatrix) "_scale_columns[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

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
      BASE aij = A->data[i];
      BASE bij = FUNCTION (gsl_spmatrix, get) (B, A->i[i], A->p[i]);
      BASE xi = FUNCTION (gsl_vector, get) (x, A->i[i]);

      if (bij != (ATOMIC) (xi * aij))
        status = 1;
    }

  gsl_test (status, NAME (gsl_spmatrix) "_scale_rows[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

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

          if (aij + bij != cij)
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

          if (aij + aij != cij)
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

          if (aij + bij != cij)
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

          if ((BASE) (aij - bij) != cij)
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

          if (bij != dij)
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
FUNCTION (test, minmax) (const size_t M, const size_t N, const int sptype,
                         const double density, gsl_rng * r)
{
  TYPE (gsl_spmatrix) * A = FUNCTION (test, random) (M, N, density, 5.0, 20.0, r);
  TYPE (gsl_spmatrix) * B;
  BASE min, max;
  size_t imin, jmin;

  FUNCTION (gsl_spmatrix, set) (A, 4, 3, (BASE) 12);
  FUNCTION (gsl_spmatrix, set) (A, 3, 5, (BASE) 5);
  FUNCTION (gsl_spmatrix, set) (A, 1, 1, (BASE) 10);
  FUNCTION (gsl_spmatrix, set) (A, 7, 3, (BASE) 2);
  FUNCTION (gsl_spmatrix, set) (A, 1, 0, (BASE) 30);

  B = FUNCTION (gsl_spmatrix, compress) (A, sptype);

  FUNCTION (gsl_spmatrix, minmax) (B, &min, &max);

  status = min != (BASE) 2;
  gsl_test (status, NAME (gsl_spmatrix) "_minmax[%zu,%zu](%s) minimum",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  status = max != (BASE) 30;
  gsl_test (status, NAME (gsl_spmatrix) "_minmax[%zu,%zu](%s) maximum",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  FUNCTION (gsl_spmatrix, min_index) (B, &imin, &jmin);

  status = imin != 7;
  gsl_test (status, NAME (gsl_spmatrix) "_min_index[%zu,%zu](%s) imin",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  status = jmin != 3;
  gsl_test (status, NAME (gsl_spmatrix) "_min_index[%zu,%zu](%s) jmin",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
}

#if !defined(UNSIGNED) && !defined(BASE_CHAR)

static void
FUNCTION (test, norm) (const size_t M, const size_t N, const int sptype)
{
  TYPE (gsl_spmatrix) * A = FUNCTION (gsl_spmatrix, alloc) (M, N);
  TYPE (gsl_spmatrix) * B;
  ATOMIC norm, norm_expected;
  size_t i, j, k = 0;

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          k++;
          FUNCTION (gsl_spmatrix, set) (A, i, j, (BASE) k);
        }
    }

  B = FUNCTION (gsl_spmatrix, compress) (A, sptype);

  norm = FUNCTION (gsl_spmatrix, norm1) (B);
  norm_expected = N*M*(M+1)/2;

  status = norm != norm_expected;
  gsl_test (status, NAME (gsl_spmatrix) "_norm1[%zu,%zu](%s)",
            M, N, FUNCTION (gsl_spmatrix, type) (B));

  FUNCTION (gsl_spmatrix, free) (A);
  FUNCTION (gsl_spmatrix, free) (B);
}

#endif

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

  FUNCTION (test, minmax) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, minmax) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, minmax) (M, N, GSL_SPMATRIX_CSR, density, r);

#if !defined(UNSIGNED) && !defined(BASE_CHAR)
  FUNCTION (test, norm) (M, N, GSL_SPMATRIX_COO);
  FUNCTION (test, norm) (M, N, GSL_SPMATRIX_CSC);
  FUNCTION (test, norm) (M, N, GSL_SPMATRIX_CSR);
#endif

  FUNCTION (test, io_ascii) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, io_ascii) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, io_ascii) (M, N, GSL_SPMATRIX_CSR, density, r);

  FUNCTION (test, io_binary) (M, N, GSL_SPMATRIX_COO, density, r);
  FUNCTION (test, io_binary) (M, N, GSL_SPMATRIX_CSC, density, r);
  FUNCTION (test, io_binary) (M, N, GSL_SPMATRIX_CSR, density, r);
}
