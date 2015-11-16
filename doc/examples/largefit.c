#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_blas.h>

/* function to be fitted */
double
func(const double t)
{
  return exp(sin(4.0 * (t - 10.0)));
}

/* construct a row of the least squares matrix */
int
build_row(const double t, gsl_vector *row)
{
  const size_t p = row->size;
  double Xj = 1.0;
  size_t j;

  for (j = 0; j < p; ++j)
    {
      gsl_vector_set(row, j, Xj);
      Xj *= t;
    }

  return 0;
}

int
solve_system(const gsl_multilarge_linear_type * T,
             const double lambda, const size_t n, const size_t p,
             gsl_vector * c)
{
  const size_t nblock = 5;         /* number of blocks to accumulate */
  const size_t nrows = n / nblock; /* number of rows per block */
  gsl_multilarge_linear_workspace * w =
    gsl_multilarge_linear_alloc(T, p);
  gsl_matrix *X = gsl_matrix_alloc(nrows, p);
  gsl_vector *y = gsl_vector_alloc(nrows);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  size_t rowidx = 0;
  double rnorm, snorm, rcond;
  double t = 10.0;
  double dt = 1.0 / (n - 1.0);

  while (rowidx < n)
    {
      size_t nleft = n - rowidx;         /* number of rows left to accumulate */
      size_t nr = GSL_MIN(nrows, nleft); /* number of rows in this block */
      gsl_matrix_view Xv = gsl_matrix_submatrix(X, 0, 0, nr, p);
      gsl_vector_view yv = gsl_vector_subvector(y, 0, nr);
      size_t i;

      /* build (X,y) block with 'nr' rows */
      for (i = 0; i < nr; ++i)
        {
          gsl_vector_view row = gsl_matrix_row(&Xv.matrix, i);
          double yi = func(t);
          double ei = gsl_ran_gaussian (r, 0.3 * yi); /* noise */

          /* construct this row of LS matrix */
          build_row(t, &row.vector);

          /* set right hand side value with added noise */
          gsl_vector_set(&yv.vector, i, yi + ei);

          t += dt;
        }

      /* accumulate (X,y) block into LS system */
      gsl_multilarge_linear_accumulate(&Xv.matrix, &yv.vector, w);

      rowidx += nr;
    }

  /* solve large LS system and store solution in c */
  gsl_multilarge_linear_solve(lambda, c, &rnorm, &snorm, w);

  /* compute reciprocal condition number */
  gsl_multilarge_linear_rcond(&rcond, w);

  fprintf(stderr, "=== Method %s ===\n", gsl_multilarge_linear_name(w));
  if (rcond != 0.0)
    fprintf(stderr, "matrix condition number = %e\n", 1.0 / rcond);
  fprintf(stderr, "residual norm  = %e\n", rnorm);
  fprintf(stderr, "solution norm  = %e\n", snorm);

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_multilarge_linear_free(w);
  gsl_rng_free(r);

  return 0;
}

int
main(int argc, char *argv[])
{
  const size_t n = 50000;   /* number of observations */
  const size_t p = 15;      /* polynomial order + 1 */
  double lambda = 0.0;      /* regularization parameter */
  gsl_vector *c_tsqr = gsl_vector_alloc(p);
  gsl_vector *c_normal = gsl_vector_alloc(p);

  if (argc > 1)
    lambda = atof(argv[1]);

  /* solve system with TSQR method */
  solve_system(gsl_multilarge_linear_tsqr, lambda, n, p, c_tsqr);

  /* solve system with Normal equations method */
  solve_system(gsl_multilarge_linear_normal, lambda, n, p, c_normal);

  /* output solutions */
  {
    gsl_vector *v = gsl_vector_alloc(p);
    double t;

    for (t = 10.0; t <= 11.0; t += 0.01)
      {
        double f_exact = func(t);
        double f_tsqr, f_normal;

        build_row(t, v);
        gsl_blas_ddot(v, c_tsqr, &f_tsqr);
        gsl_blas_ddot(v, c_normal, &f_normal);

        printf("%f %e %e %e\n", t, f_exact, f_tsqr, f_normal);
      }

    gsl_vector_free(v);
  }

  gsl_vector_free(c_tsqr);
  gsl_vector_free(c_normal);

  return 0;
}
