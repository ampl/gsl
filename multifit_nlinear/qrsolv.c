/* This function computes the solution to the least squares system

   phi = [ A x =  b , lambda D x = 0 ]^2
    
   where A is an M by N matrix, D is an N by N diagonal matrix, lambda
   is a scalar parameter and b is a vector of length M.

   The function requires the factorization of A into A = Q R P^T,
   where Q is an orthogonal matrix, R is an upper triangular matrix
   with diagonal elements of non-increasing magnitude and P is a
   permuation matrix. The system above is then equivalent to

   [ R z = Q^T b, P^T (lambda D) P z = 0 ]

   where x = P z. If this system does not have full rank then a least
   squares solution is obtained.  On output the function also provides
   an upper triangular matrix S such that

   P^T (A^T A + lambda^2 D^T D) P = S^T S

   Parameters,
   
   r: On input, contains the full upper triangle of R. The diagonal
   elements are modified but restored on output. The full
   upper triangle of R is not modified.

   p: the encoded form of the permutation matrix P. column j of P is
   column p[j] of the identity matrix.

   lambda, diag: contains the scalar lambda and the diagonal elements
   of the matrix D

   qtb: contains the product Q^T b

   S: on output contains the matrix S, n-by-n

   x: on output contains the least squares solution of the system

   work: is a workspace of length N

   */

static int
qrsolv (gsl_matrix * r, const gsl_permutation * p, const double lambda,
        const gsl_vector * diag, const gsl_vector * qtb,
        gsl_matrix * S, gsl_vector * x, gsl_vector * work)
{
  size_t n = r->size2;

  size_t i, j, k, nsing;

  /* Copy r and qtb to preserve input and initialise s. In particular,
     save the diagonal elements of r in x */

  for (j = 0; j < n; j++)
    {
      double rjj = gsl_matrix_get (r, j, j);
      double qtbj = gsl_vector_get (qtb, j);

      for (i = j + 1; i < n; i++)
        {
          double rji = gsl_matrix_get (r, j, i);
          gsl_matrix_set (S, i, j, rji);
        }

      gsl_vector_set (x, j, rjj);
      gsl_vector_set (work, j, qtbj);
    }

  /* Eliminate the diagonal matrix d using a Givens rotation */

  for (j = 0; j < n; j++)
    {
      double qtbpj;

      size_t pj = gsl_permutation_get (p, j);

      double diagpj = lambda * gsl_vector_get (diag, pj);

      if (diagpj == 0)
        {
          continue;
        }

      gsl_matrix_set (S, j, j, diagpj);

      for (k = j + 1; k < n; k++)
        {
          gsl_matrix_set (S, k, k, 0.0);
        }

      /* The transformations to eliminate the row of d modify only a
         single element of qtb beyond the first n, which is initially
         zero */

      qtbpj = 0;

      for (k = j; k < n; k++)
        {
          /* Determine a Givens rotation which eliminates the
             appropriate element in the current row of d */

          double sine, cosine;

          double wk = gsl_vector_get (work, k);
          double rkk = gsl_matrix_get (r, k, k);
          double skk = gsl_matrix_get (S, k, k);

          if (skk == 0)
            {
              continue;
            }

          if (fabs (rkk) < fabs (skk))
            {
              double cotangent = rkk / skk;
              sine = 0.5 / sqrt (0.25 + 0.25 * cotangent * cotangent);
              cosine = sine * cotangent;
            }
          else
            {
              double tangent = skk / rkk;
              cosine = 0.5 / sqrt (0.25 + 0.25 * tangent * tangent);
              sine = cosine * tangent;
            }

          /* Compute the modified diagonal element of r and the
             modified element of [qtb,0] */

          {
            double new_rkk = cosine * rkk + sine * skk;
            double new_wk = cosine * wk + sine * qtbpj;
            
            qtbpj = -sine * wk + cosine * qtbpj;

            gsl_matrix_set(r, k, k, new_rkk);
            gsl_matrix_set(S, k, k, new_rkk);
            gsl_vector_set(work, k, new_wk);
          }

          /* Accumulate the transformation in the row of s */

          for (i = k + 1; i < n; i++)
            {
              double sik = gsl_matrix_get (S, i, k);
              double sii = gsl_matrix_get (S, i, i);
              
              double new_sik = cosine * sik + sine * sii;
              double new_sii = -sine * sik + cosine * sii;

              gsl_matrix_set(S, i, k, new_sik);
              gsl_matrix_set(S, i, i, new_sii);
            }
        }

      /* Store the corresponding diagonal element of s and restore the
         corresponding diagonal element of r */

      {
        double xj = gsl_vector_get(x, j);
        gsl_matrix_set (r, j, j, xj);
      }

    }

  /* Solve the triangular system for z. If the system is singular then
     obtain a least squares solution */

  nsing = n;

  for (j = 0; j < n; j++)
    {
      double sjj = gsl_matrix_get (S, j, j);

      if (sjj == 0)
        {
          nsing = j;
          break;
        }
    }

  for (j = nsing; j < n; j++)
    {
      gsl_vector_set (work, j, 0.0);
    }

  for (k = 0; k < nsing; k++)
    {
      double sum = 0;

      j = (nsing - 1) - k;

      for (i = j + 1; i < nsing; i++)
        {
          sum += gsl_matrix_get(S, i, j) * gsl_vector_get(work, i);
        }

      {
        double wj = gsl_vector_get (work, j);
        double sjj = gsl_matrix_get (S, j, j);

        gsl_vector_set (work, j, (wj - sum) / sjj);
      }
    }

  /* Permute the components of z back to the components of x */

  for (j = 0; j < n; j++)
    {
      size_t pj = gsl_permutation_get (p, j);
      double wj = gsl_vector_get (work, j);

      gsl_vector_set (x, pj, wj);
    }

  return GSL_SUCCESS;
}
