/* linalg/cholesky_common.c
 * 
 * Copyright (C) 2016 Patrick Alken
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

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

/*
cholesky_swap_rowcol()
  Swap rows and columns i and j of symmetric matrix A, updating only the
lower triangle

Inputs: A - symmetric matrix A, stored in lower triangle
        i - index i
        j - index j

Return: success/error
*/

static int
cholesky_swap_rowcol(gsl_matrix * A, const size_t i, const size_t j)
{
  if (i != j)
    {
      const size_t N = A->size1;
      double *Aii, *Ajj;
      size_t ii, jj, k;

      /* organize so that ii < jj */
      if (i < j)
        {
          ii = i;
          jj = j;
        }
      else
        {
          ii = j;
          jj = i;
        }

      /* swap subrows A(i,1:i-1) with A(j,1:i-1) */
      for (k = 0; k < ii; ++k)
        {
          double *Aik = gsl_matrix_ptr(A, ii, k);
          double *Ajk = gsl_matrix_ptr(A, jj, k);
          SWAP(*Aik, *Ajk);
        }

      /* swap subrow A(j,i+1:j-1) with subcolumn A(i+1:j-1,i) */
      for (k = ii + 1; k < jj; ++k)
        {
          double *Ajk = gsl_matrix_ptr(A, jj, k);
          double *Aki = gsl_matrix_ptr(A, k, ii);
          SWAP(*Ajk, *Aki);
        }

      /* swap subcolumns A(j+1:N,i) with A(j+1:N,j) */
      for (k = jj + 1; k < N; ++k)
        {
          double *Aki = gsl_matrix_ptr(A, k, ii);
          double *Akj = gsl_matrix_ptr(A, k, jj);
          SWAP(*Aki, *Akj);
        }

      /* now swap diagonal elements A(i,i) and A(j,j) */
      Aii = gsl_matrix_ptr(A, ii, ii);
      Ajj = gsl_matrix_ptr(A, jj, jj);
      SWAP(*Aii, *Ajj);
    }

  return GSL_SUCCESS;
}
