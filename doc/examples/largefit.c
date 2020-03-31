#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

/* function to be fitted */
double
func(const double t)
{
  double x = sin(10.0 * t);
  return exp(x*x*x);
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
solve_system(const int print_data, const gsl_multilarge_linear_type * T,
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
  const size_t nlcurve = 200;
  gsl_vector *reg_param = gsl_vector_alloc(nlcurve);
  gsl_vector *rho = gsl_vector_calloc(nlcurve);
  gsl_vector *eta = gsl_vector_calloc(nlcurve);
  size_t rowidx = 0;
  double rnorm, snorm, rcond;
  double t = 0.0;
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
          double fi = func(t);
          double ei = gsl_ran_gaussian (r, 0.1 * fi); /* noise */
          double yi = fi + ei;

          /* construct this row of LS matrix */
          build_row(t, &row.vector);

          /* set right hand side value with added noise */
          gsl_vector_set(&yv.vector, i, yi);

          if (print_data && (i % 100 == 0))
            printf("%f %f\n", t, yi);

          t += dt;
        }

      /* accumulate (X,y) block into LS system */
      gsl_multilarge_linear_accumulate(&Xv.matrix, &yv.vector, w);

      rowidx += nr;
    }

  if (print_data)
    printf("\n\n");

  /* compute L-curve */
  gsl_multilarge_linear_lcurve(reg_param, rho, eta, w);

  /* solve large LS system and store solution in c */
  gsl_multilarge_linear_solve(lambda, c, &rnorm, &snorm, w);

  /* compute reciprocal condition number */
  gsl_multilarge_linear_rcond(&rcond, w);

  fprintf(stderr, "=== Method %s ===\n", gsl_multilarge_linear_name(w));
  fprintf(stderr, "condition number = %e\n", 1.0 / rcond);
  fprintf(stderr, "residual norm    = %e\n", rnorm);
  fprintf(stderr, "solution norm    = %e\n", snorm);

  /* output L-curve */
  {
    size_t i;
    for (i = 0; i < nlcurve; ++i)
      {
        printf("%.12e %.12e %.12e\n",
               gsl_vector_get(reg_param, i),
               gsl_vector_get(rho, i),
               gsl_vector_get(eta, i));
      }
    printf("\n\n");
  }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_multilarge_linear_free(w);
  gsl_rng_free(r);
  gsl_vector_free(reg_param);
  gsl_vector_free(rho);
  gsl_vector_free(eta);

  return 0;
}

int
main(int argc, char *argv[])
{
  const size_t n = 50000;   /* number of observations */
  const size_t p = 16;      /* polynomial order + 1 */
  double lambda = 0.0;      /* regularization parameter */
  gsl_vector *c_tsqr = gsl_vector_calloc(p);
  gsl_vector *c_normal = gsl_vector_calloc(p);

  if (argc > 1)
    lambda = atof(argv[1]);

  /* turn off error handler so normal equations method won't abort */
  gsl_set_error_handler_off();

  /* solve system with TSQR method */
  solve_system(1, gsl_multilarge_linear_tsqr, lambda, n, p, c_tsqr);

  /* solve system with Normal equations method */
  solve_system(0, gsl_multilarge_linear_normal, lambda, n, p, c_normal);

  /* output solutions */
  {
    gsl_vector *v = gsl_vector_alloc(p);
    double t;

    for (t = 0.0; t <= 1.0; t += 0.01)
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
