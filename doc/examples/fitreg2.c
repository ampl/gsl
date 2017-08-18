#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

static int
hilbert_matrix(gsl_matrix * m)
{
  const size_t N = m->size1;
  const size_t M = m->size2;
  size_t i, j;

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < M; j++)
        {
          gsl_matrix_set(m, i, j, 1.0/(i+j+1.0));
        }
    }

  return GSL_SUCCESS;
}

int
main()
{
  const size_t n = 10; /* number of observations */
  const size_t p = 8;  /* number of model parameters */
  size_t i;
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *y = gsl_vector_alloc(n);

  /* construct Hilbert matrix and rhs vector */
  hilbert_matrix(X);

  {
    double val = 1.0;
    for (i = 0; i < n; ++i)
      {
        gsl_vector_set(y, i, val);
        val *= -1.0;
      }
  }

  {
    const size_t npoints = 200;                   /* number of points on L-curve and GCV curve */
    gsl_multifit_linear_workspace *w =
      gsl_multifit_linear_alloc(n, p);
    gsl_vector *c = gsl_vector_alloc(p);          /* OLS solution */
    gsl_vector *c_lcurve = gsl_vector_alloc(p);   /* regularized solution (L-curve) */
    gsl_vector *c_gcv = gsl_vector_alloc(p);      /* regularized solution (GCV) */
    gsl_vector *reg_param = gsl_vector_alloc(npoints);
    gsl_vector *rho = gsl_vector_alloc(npoints);  /* residual norms */
    gsl_vector *eta = gsl_vector_alloc(npoints);  /* solution norms */
    gsl_vector *G = gsl_vector_alloc(npoints);    /* GCV function values */
    double lambda_l;                              /* optimal regularization parameter (L-curve) */
    double lambda_gcv;                            /* optimal regularization parameter (GCV) */
    double G_gcv;                                 /* G(lambda_gcv) */
    size_t reg_idx;                               /* index of optimal lambda */
    double rcond;                                 /* reciprocal condition number of X */
    double chisq, rnorm, snorm;

    /* compute SVD of X */
    gsl_multifit_linear_svd(X, w);

    rcond = gsl_multifit_linear_rcond(w);
    fprintf(stderr, "matrix condition number = %e\n", 1.0 / rcond);

    /* unregularized (standard) least squares fit, lambda = 0 */
    gsl_multifit_linear_solve(0.0, X, y, c, &rnorm, &snorm, w);
    chisq = pow(rnorm, 2.0);

    fprintf(stderr, "\n=== Unregularized fit ===\n");
    fprintf(stderr, "residual norm = %g\n", rnorm);
    fprintf(stderr, "solution norm = %g\n", snorm);
    fprintf(stderr, "chisq/dof = %g\n", chisq / (n - p));

    /* calculate L-curve and find its corner */
    gsl_multifit_linear_lcurve(y, reg_param, rho, eta, w);
    gsl_multifit_linear_lcorner(rho, eta, &reg_idx);

    /* store optimal regularization parameter */
    lambda_l = gsl_vector_get(reg_param, reg_idx);

    /* regularize with lambda_l */
    gsl_multifit_linear_solve(lambda_l, X, y, c_lcurve, &rnorm, &snorm, w);
    chisq = pow(rnorm, 2.0) + pow(lambda_l * snorm, 2.0);

    fprintf(stderr, "\n=== Regularized fit (L-curve) ===\n");
    fprintf(stderr, "optimal lambda: %g\n", lambda_l);
    fprintf(stderr, "residual norm = %g\n", rnorm);
    fprintf(stderr, "solution norm = %g\n", snorm);
    fprintf(stderr, "chisq/dof = %g\n", chisq / (n - p));

    /* calculate GCV curve and find its minimum */
    gsl_multifit_linear_gcv(y, reg_param, G, &lambda_gcv, &G_gcv, w);

    /* regularize with lambda_gcv */
    gsl_multifit_linear_solve(lambda_gcv, X, y, c_gcv, &rnorm, &snorm, w);
    chisq = pow(rnorm, 2.0) + pow(lambda_gcv * snorm, 2.0);

    fprintf(stderr, "\n=== Regularized fit (GCV) ===\n");
    fprintf(stderr, "optimal lambda: %g\n", lambda_gcv);
    fprintf(stderr, "residual norm = %g\n", rnorm);
    fprintf(stderr, "solution norm = %g\n", snorm);
    fprintf(stderr, "chisq/dof = %g\n", chisq / (n - p));

    /* output L-curve and GCV curve */
    for (i = 0; i < npoints; ++i)
      {
        printf("%e %e %e %e\n",
               gsl_vector_get(reg_param, i),
               gsl_vector_get(rho, i),
               gsl_vector_get(eta, i),
               gsl_vector_get(G, i));
      }

    /* output L-curve corner point */
    printf("\n\n%f %f\n",
           gsl_vector_get(rho, reg_idx),
           gsl_vector_get(eta, reg_idx));

    /* output GCV curve corner minimum */
    printf("\n\n%e %e\n",
           lambda_gcv,
           G_gcv);

    gsl_multifit_linear_free(w);
    gsl_vector_free(c);
    gsl_vector_free(c_lcurve);
    gsl_vector_free(reg_param);
    gsl_vector_free(rho);
    gsl_vector_free(eta);
    gsl_vector_free(G);
  }

  gsl_matrix_free(X);
  gsl_vector_free(y);

  return 0;
}
