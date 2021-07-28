/* spmatrix/oper_complex_source.c
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

static size_t
FUNCTION (spmatrix, scatter) (const TYPE (gsl_spmatrix) * A, const size_t j, int * w,
                              ATOMIC * x, const int mark, TYPE (gsl_spmatrix) * C, size_t nz);

int
FUNCTION (gsl_spmatrix, scale) (TYPE (gsl_spmatrix) * m, const BASE x)
{
  const ATOMIC xr = GSL_REAL (x);
  const ATOMIC xi = GSL_IMAG (x);
  size_t i;

  for (i = 0; i < m->nz; ++i)
    {
      ATOMIC mr = m->data[2 * i];
      ATOMIC mi = m->data[2 * i + 1];

      m->data[2 * i] = mr * xr - mi * xi;
      m->data[2 * i + 1] = mi * xr + mr * xi;
    }

  return GSL_SUCCESS;
}

/* m := m * diag(x) */
int
FUNCTION (gsl_spmatrix, scale_columns) (TYPE (gsl_spmatrix) * m, const TYPE (gsl_vector) * x)
{
  if (m->size2 != x->size)
    {
      GSL_ERROR("x vector length does not match matrix", GSL_EBADLEN);
    }
  else
    {
      ATOMIC * Ad = m->data;

      if (GSL_SPMATRIX_ISCSC(m))
        {
          const int *Ap = m->p;
          int p;
          size_t j;

          for (j = 0; j < m->size2; ++j)
            {
              BASE xj = FUNCTION (gsl_vector, get) (x, j);
              ATOMIC xr = GSL_REAL(xj);
              ATOMIC xi = GSL_IMAG(xj);

              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  ATOMIC ar = Ad[2 * p];
                  ATOMIC ai = Ad[2 * p + 1];

                  Ad[2 * p] = ar * xr - ai * xi;
                  Ad[2 * p + 1] = ai * xr + ar * xi;
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSR(m))
        {
          const int *Aj = m->i;
          size_t i;

          for (i = 0; i < m->nz; ++i)
            {
              BASE y = FUNCTION (gsl_vector, get) (x, Aj[i]);
              ATOMIC yr = GSL_REAL(y);
              ATOMIC yi = GSL_IMAG(y);
              ATOMIC ar = Ad[2 * i];
              ATOMIC ai = Ad[2 * i + 1];

              Ad[2 * i] = ar * yr - ai * yi;
              Ad[2 * i + 1] = ai * yr + ar * yi;
            }
        }
      else if (GSL_SPMATRIX_ISCOO(m))
        {
          const int *Aj = m->p;
          size_t i;

          for (i = 0; i < m->nz; ++i)
            {
              BASE y = FUNCTION (gsl_vector, get) (x, Aj[i]);
              ATOMIC yr = GSL_REAL(y);
              ATOMIC yi = GSL_IMAG(y);
              ATOMIC ar = Ad[2 * i];
              ATOMIC ai = Ad[2 * i + 1];

              Ad[2 * i] = ar * yr - ai * yi;
              Ad[2 * i + 1] = ai * yr + ar * yi;
            }
        }
      else
        {
          GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
        }

      return GSL_SUCCESS;
    }
}

/* m := diag(x) * m */
int
FUNCTION (gsl_spmatrix, scale_rows) (TYPE (gsl_spmatrix) * m, const TYPE (gsl_vector) * x)
{
  if (m->size1 != x->size)
    {
      GSL_ERROR("x vector length does not match matrix", GSL_EBADLEN);
    }
  else
    {
      ATOMIC * Ad = m->data;

      if (GSL_SPMATRIX_ISCSC(m))
        {
          const int *Ai = m->i;
          size_t i;

          for (i = 0; i < m->nz; ++i)
            {
              BASE y = FUNCTION (gsl_vector, get) (x, Ai[i]);
              ATOMIC yr = GSL_REAL(y);
              ATOMIC yi = GSL_IMAG(y);
              ATOMIC ar = Ad[2 * i];
              ATOMIC ai = Ad[2 * i + 1];

              Ad[2 * i] = ar * yr - ai * yi;
              Ad[2 * i + 1] = ai * yr + ar * yi;
            }
        }
      else if (GSL_SPMATRIX_ISCSR(m))
        {
          const int *Ap = m->p;
          int p;
          size_t i;

          for (i = 0; i < m->size1; ++i)
            {
              BASE y = FUNCTION (gsl_vector, get) (x, i);
              ATOMIC yr = GSL_REAL(y);
              ATOMIC yi = GSL_IMAG(y);

              for (p = Ap[i]; p < Ap[i + 1]; ++p)
                {
                  ATOMIC ar = Ad[2 * p];
                  ATOMIC ai = Ad[2 * p + 1];
                  Ad[2 * p] = ar * yr - ai * yi;
                  Ad[2 * p + 1] = ai * yr + ar * yi;
                }
            }
        }
      else if (GSL_SPMATRIX_ISCOO(m))
        {
          const int *Ai = m->i;
          size_t i;

          for (i = 0; i < m->nz; ++i)
            {
              BASE y = FUNCTION (gsl_vector, get) (x, Ai[i]);
              ATOMIC yr = GSL_REAL(y);
              ATOMIC yi = GSL_IMAG(y);
              ATOMIC ar = Ad[2 * i];
              ATOMIC ai = Ad[2 * i + 1];

              Ad[2 * i] = ar * yr - ai * yi;
              Ad[2 * i + 1] = ai * yr + ar * yi;
            }
        }
      else
        {
          GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_spmatrix_add()
  Add two sparse matrices

Inputs: c - (output) a + b
        a - (input) sparse matrix
        b - (input) sparse matrix

Return: success or error
*/

int
FUNCTION (gsl_spmatrix, add) (TYPE (gsl_spmatrix) * c, const TYPE (gsl_spmatrix) * a,
                              const TYPE (gsl_spmatrix) * b)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (b->size1 != M || b->size2 != N || c->size1 != M || c->size2 != N)
    {
      GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
    }
  else if (a->sptype != b->sptype || a->sptype != c->sptype)
    {
      GSL_ERROR("matrices must have same sparse storage format",
                GSL_EINVAL);
    }
  else if (GSL_SPMATRIX_ISCOO(a))
    {
      GSL_ERROR("COO format not yet supported", GSL_EINVAL);
    }
  else
    {
      int status = GSL_SUCCESS;
      int * w = a->work.work_int;
      ATOMIC * x = c->work.work_atomic;
      int *Cp, *Ci;
      ATOMIC * Cd;
      int p;
      size_t j;
      size_t nz = 0; /* number of non-zeros in c */
      size_t inner_size, outer_size;

      if (GSL_SPMATRIX_ISCSC(a))
        {
          inner_size = M;
          outer_size = N;
        }
      else if (GSL_SPMATRIX_ISCSR(a))
        {
          inner_size = N;
          outer_size = M;
        }
      else
        {
          GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
        }

      if (c->nzmax < a->nz + b->nz)
        {
          status = FUNCTION (gsl_spmatrix, realloc) (a->nz + b->nz, c);
          if (status)
            return status;
        }

      /* initialize w = 0 */
      for (j = 0; j < inner_size; ++j)
        w[j] = 0;

      Ci = c->i;
      Cp = c->p;
      Cd = c->data;

      for (j = 0; j < outer_size; ++j)
        {
          Cp[j] = nz;

          /* CSC: x += A(:,j); CSR: x += A(j,:) */
          nz = FUNCTION (spmatrix, scatter) (a, j, w, x, (int) (j + 1), c, nz);

          /* CSC: x += B(:,j); CSR: x += B(j,:) */
          nz = FUNCTION (spmatrix, scatter) (b, j, w, x, (int) (j + 1), c, nz);

          for (p = Cp[j]; p < (int) nz; ++p)
            {
              Cd[2 * p] = x[2 * Ci[p]];
              Cd[2 * p + 1] = x[2 * Ci[p] + 1];
            }
        }

      /* finalize last column of c */
      Cp[j] = nz;
      c->nz = nz;

      return status;
    }
}

/*
gsl_spmatrix_dense_add()
  a := a + b

Inputs: a - (input/output) on input, dense matrix a; on output, a + b
        b - (input) sparse matrix

Return: success or error
*/

int
FUNCTION (gsl_spmatrix, dense_add) (TYPE (gsl_matrix) * a, const TYPE (gsl_spmatrix) * b)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (b->size1 != M || b->size2 != N)
    {
      GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
    }
  else
    {
      const size_t tda_a = a->tda;
      const ATOMIC * bd = b->data;

      /* check for quick return */
      if (b->nz == 0)
        return GSL_SUCCESS;

      if (GSL_SPMATRIX_ISCOO(b))
        {
          const int * bi = b->i;
          const int * bj = b->p;
          size_t n;

          for (n = 0; n < b->nz; ++n)
            {
              const size_t idx = 2 * (bi[n] * tda_a + bj[n]);
              a->data[idx] += bd[2 * n];
              a->data[idx + 1] += bd[2 * n + 1];
            }
        }
      else if (GSL_SPMATRIX_ISCSC(b))
        {
          const int * bi = b->i;
          const int * bp = b->p;
          size_t j;
          int p;

          for (j = 0; j < N; ++j)
            {
              for (p = bp[j]; p < bp[j + 1]; ++p)
                {
                  const size_t idx = 2 * (bi[p] * tda_a + j);
                  a->data[idx] += bd[2 * p];
                  a->data[idx + 1] += bd[2 * p + 1];
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSR(b))
        {
          const int * bj = b->i;
          const int * bp = b->p;
          size_t i;
          int p;

          for (i = 0; i < M; ++i)
            {
              for (p = bp[i]; p < bp[i + 1]; ++p)
                {
                  const size_t idx = 2 * (i * tda_a + bj[p]);
                  a->data[idx] += bd[2 * p];
                  a->data[idx + 1] += bd[2 * p + 1];
                }
            }
        }

      return GSL_SUCCESS;
    }
}

#ifndef GSL_DISABLE_DEPRECATED

int
FUNCTION (gsl_spmatrix, add_to_dense) (TYPE (gsl_matrix) * a, const TYPE (gsl_spmatrix) * b)
{
  return FUNCTION (gsl_spmatrix, dense_add) (a, b);
}

#endif

/*
gsl_spmatrix_dense_sub()
  a := a - b

Inputs: a - (input/output) on input, dense matrix a; on output, a - b
        b - (input) sparse matrix

Return: success or error
*/

int
FUNCTION (gsl_spmatrix, dense_sub) (TYPE (gsl_matrix) * a, const TYPE (gsl_spmatrix) * b)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (b->size1 != M || b->size2 != N)
    {
      GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
    }
  else
    {
      const size_t tda_a = a->tda;
      const ATOMIC * bd = b->data;

      /* check for quick return */
      if (b->nz == 0)
        return GSL_SUCCESS;

      if (GSL_SPMATRIX_ISCOO(b))
        {
          const int * bi = b->i;
          const int * bj = b->p;
          size_t n;

          for (n = 0; n < b->nz; ++n)
            {
              const size_t idx = 2 * (bi[n] * tda_a + bj[n]);
              a->data[idx] -= bd[2 * n];
              a->data[idx + 1] -= bd[2 * n + 1];
            }
        }
      else if (GSL_SPMATRIX_ISCSC(b))
        {
          const int * bi = b->i;
          const int * bp = b->p;
          size_t j;
          int p;

          for (j = 0; j < N; ++j)
            {
              for (p = bp[j]; p < bp[j + 1]; ++p)
                {
                  const size_t idx = 2 * (bi[p] * tda_a + j);
                  a->data[idx] -= bd[2 * p];
                  a->data[idx + 1] -= bd[2 * p + 1];
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSR(b))
        {
          const int * bj = b->i;
          const int * bp = b->p;
          size_t i;
          int p;

          for (i = 0; i < M; ++i)
            {
              for (p = bp[i]; p < bp[i + 1]; ++p)
                {
                  const size_t idx = 2 * (i * tda_a + bj[p]);
                  a->data[idx] -= bd[2 * p];
                  a->data[idx + 1] -= bd[2 * p + 1];
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_spmatrix_d2sp()
  Convert a dense gsl_matrix to sparse (COO) format

Inputs: T - (output) sparse matrix in COO format
        A - (input) dense matrix to convert
*/

int
FUNCTION (gsl_spmatrix, d2sp) (TYPE (gsl_spmatrix) * T, const TYPE (gsl_matrix) * A)
{
  if (T->size1 != A->size1 || T->size2 != A->size2)
    {
      GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
    }
  else if (!GSL_SPMATRIX_ISCOO(T))
    {
      GSL_ERROR ("sparse matrix must be in COO format", GSL_EINVAL);
    }
  else
    {
      size_t i, j;

      FUNCTION (gsl_spmatrix, set_zero) (T);

      for (i = 0; i < A->size1; ++i)
        {
          for (j = 0; j < A->size2; ++j)
            {
              BASE x = FUNCTION (gsl_matrix, get) (A, i, j);

              if (!(GSL_REAL(x) == (ATOMIC) 0 && GSL_IMAG(x) == (ATOMIC) 0))
                FUNCTION (gsl_spmatrix, set) (T, i, j, x);
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_spmatrix_sp2d()
  Convert a sparse matrix to dense format
*/

int
FUNCTION (gsl_spmatrix, sp2d) (TYPE (gsl_matrix) * A, const TYPE (gsl_spmatrix) * S)
{
  if (A->size1 != S->size1 || A->size2 != S->size2)
    {
      GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
    }
  else
    {
      FUNCTION (gsl_matrix, set_zero) (A);

      if (GSL_SPMATRIX_ISCOO(S))
        {
          size_t n;

          for (n = 0; n < S->nz; ++n)
            {
              int i = S->i[n];
              int j = S->p[n];
              BASE x = *(BASE *) &S->data[2 * n];

              FUNCTION (gsl_matrix, set) (A, i, j, x);
            }
        }
      else if (GSL_SPMATRIX_ISCSR(S))
        {
          const int *Sj = S->i;
          const int *Sp = S->p;
          const ATOMIC *Sd = S->data;
          size_t i;
          int p;

          for (i = 0; i < S->size1; ++i)
            {
              for (p = Sp[i]; p < Sp[i + 1]; ++p)
                {
                  BASE x;
                  GSL_SET_COMPLEX(&x, Sd[2*p], Sd[2*p + 1]);
                  FUNCTION (gsl_matrix, set) (A, i, Sj[p], x);
                }
            }
        }
      else if (GSL_SPMATRIX_ISCSC(S))
        {
          const int *Si = S->i;
          const int *Sp = S->p;
          const ATOMIC *Sd = S->data;
          size_t j;
          int p;

          for (j = 0; j < S->size2; ++j)
            {
              for (p = Sp[j]; p < Sp[j + 1]; ++p)
                {
                  BASE x;
                  GSL_SET_COMPLEX(&x, Sd[2*p], Sd[2*p+1]);
                  FUNCTION (gsl_matrix, set) (A, Si[p], j, x);
                }
            }
        }
      else
        {
          GSL_ERROR("unknown sparse matrix type", GSL_EINVAL);
        }

      return GSL_SUCCESS;
    }
}

/*
spmatrix_scatter()

  Keep a running total x -> x + alpha*A(:,j) for adding matrices together in CSC,
which will eventually be stored in C(:,j)

  When a new non-zero element with row index i is found, update C->i with
the row index. C->data is updated only by the calling function after all
matrices have been added via this function.

Inputs: A     - sparse matrix m-by-n
        j     - column index
        w     - keeps track which rows of column j have been added to C;
                initialize to 0 prior to first call
        x     - column vector of length m
        mark  -
        C     - output matrix whose jth column will be added to A(:,j)
        nz    - (input/output) number of non-zeros in matrix C

Notes:
1) This function is designed to be called successively when adding multiple
matrices together. Column j of C is stored contiguously as per CSC but not
necessarily in order - ie: the row indices C->i may not be in ascending order.

2) based on CSparse routine cs_scatter
*/

static size_t
FUNCTION (spmatrix, scatter) (const TYPE (gsl_spmatrix) * A, const size_t j, int * w,
                              ATOMIC * x, const int mark, TYPE (gsl_spmatrix) * C, size_t nz)
{
  int p;
  int * Ai = A->i;
  int * Ap = A->p;
  ATOMIC * Ad = A->data;
  int * Ci = C->i;

  for (p = Ap[j]; p < Ap[j + 1]; ++p)
    {
      int i = Ai[p];             /* A(i,j) is nonzero */

      if (w[i] < mark)           /* check if row i has been stored in column j yet */
        {
          w[i] = mark;           /* i is new entry in column j */
          Ci[nz++] = i;          /* add i to pattern of C(:,j) */
          x[2*i] = Ad[2*p];      /* x(i) = A(i,j) */
          x[2*i+1] = Ad[2*p+1];
        }
      else                       /* this (i,j) exists in C from a previous call */
        {
          x[2*i] += Ad[2*p];     /* add A(i,j) to C(i,j) */
          x[2*i+1] += Ad[2*p+1];
        }
    }

  return (nz) ;
}
