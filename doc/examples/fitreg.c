#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>

int
main()
{
  const size_t n = 1000; /* number of observations */
  const size_t p = 2;    /* number of model parameters */
  size_t i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *y = gsl_vector_alloc(n);

  for (i = 0; i < n; ++i)
    {
      /* generate first random variable u */
      double ui = 5.0 * gsl_ran_gaussian(r, 1.0);

      /* set v = u + noise */
      double vi = ui + gsl_ran_gaussian(r, 0.001);

      /* set y = u + v + noise */
      double yi = ui + vi + gsl_ran_gaussian(r, 1.0);

      /* since u =~ v, the matrix X is ill-conditioned */
      gsl_matrix_set(X, i, 0, ui);
      gsl_matrix_set(X, i, 1, vi);

      /* rhs vector */
      gsl_vector_set(y, i, yi);
    }

  {
    const size_t nL = 200;                   /* number of points on L-curve */
    gsl_multifit_linear_workspace *w =
      gsl_multifit_linear_alloc(n, p);
    gsl_vector *c = gsl_vector_alloc(p);     /* OLS solution */
    gsl_vector *c_reg = gsl_vector_alloc(p); /* regularized solution */
    gsl_vector *reg_param = gsl_vector_alloc(nL);
    gsl_vector *rho = gsl_vector_alloc(nL);  /* residual norms */
    gsl_vector *eta = gsl_vector_alloc(nL);  /* solution norms */
    double lambda;                           /* optimal regularization parameter */
    size_t reg_idx;                          /* index of optimal lambda */
    double rcond;                            /* reciprocal condition number of X */
    double chisq, rnorm, snorm;

    /* compute SVD of X */
    gsl_multifit_linear_svd(X, w);

    rcond = gsl_multifit_linear_rcond(w);
    fprintf(stderr, "matrix condition number = %e\n", 1.0 / rcond);

    /* unregularized (standard) least squares fit, lambda = 0 */
    gsl_multifit_linear_solve(0.0, X, y, c, &rnorm, &snorm, w);
    chisq = pow(rnorm, 2.0);

    fprintf(stderr, "=== Unregularized fit ===\n");
    fprintf(stderr, "best fit: y = %g u + %g v\n",
      gsl_vector_get(c, 0), gsl_vector_get(c, 1));
    fprintf(stderr, "chisq/dof = %g\n", chisq / (n - p));

    /* calculate L-curve and find its corner */
    gsl_multifit_linear_lcurve(y, reg_param, rho, eta, w);
    gsl_multifit_linear_lcorner(rho, eta, &reg_idx);

    /* store optimal regularization parameter */
    lambda = gsl_vector_get(reg_param, reg_idx);

    /* output L-curve */
    for (i = 0; i < nL; ++i)
      printf("%f %f\n", gsl_vector_get(rho, i), gsl_vector_get(eta, i));

    /* output L-curve corner point */
    printf("\n\n%f %f\n",
           gsl_vector_get(rho, reg_idx),
           gsl_vector_get(eta, reg_idx));

    /* regularize with lambda */
    gsl_multifit_linear_solve(lambda, X, y, c_reg, &rnorm, &snorm, w);
    chisq = pow(rnorm, 2.0) + pow(lambda * snorm, 2.0);

    fprintf(stderr, "=== Regularized fit ===\n");
    fprintf(stderr, "optimal lambda: %g\n", lambda);
    fprintf(stderr, "best fit: y = %g u + %g v\n",
            gsl_vector_get(c_reg, 0), gsl_vector_get(c_reg, 1));
    fprintf(stderr, "chisq/dof = %g\n", chisq / (n - p));

    gsl_multifit_linear_free(w);
    gsl_vector_free(c);
    gsl_vector_free(c_reg);
    gsl_vector_free(reg_param);
    gsl_vector_free(rho);
    gsl_vector_free(eta);
  }

  gsl_rng_free(r);
  gsl_matrix_free(X);
  gsl_vector_free(y);

  return 0;
}
