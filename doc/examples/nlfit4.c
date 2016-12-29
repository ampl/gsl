#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_spmatrix.h>

/* parameters for functions */
struct model_params
{
  double alpha;
  gsl_spmatrix *J;
};

/* penalty function */
int
penalty_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct model_params *par = (struct model_params *) params;
  const double sqrt_alpha = sqrt(par->alpha);
  const size_t p = x->size;
  size_t i;
  double sum = 0.0;

  for (i = 0; i < p; ++i)
    {
      double xi = gsl_vector_get(x, i);

      gsl_vector_set(f, i, sqrt_alpha*(xi - 1.0));

      sum += xi * xi;
    }

  gsl_vector_set(f, p, sum - 0.25);

  return GSL_SUCCESS;
}

int
penalty_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
            const gsl_vector * u, void * params, gsl_vector * v,
            gsl_matrix * JTJ)
{
  struct model_params *par = (struct model_params *) params;
  const size_t p = x->size;
  size_t j;

  /* store 2*x in last row of J */
  for (j = 0; j < p; ++j)
    {
      double xj = gsl_vector_get(x, j);
      gsl_spmatrix_set(par->J, p, j, 2.0 * xj);
    }

  /* compute v = op(J) u */
  if (v)
    gsl_spblas_dgemv(TransJ, 1.0, par->J, u, 0.0, v);

  if (JTJ)
    {
      gsl_vector_view diag = gsl_matrix_diagonal(JTJ);

      /* compute J^T J = [ alpha*I_p + 4 x x^T ] */
      gsl_matrix_set_zero(JTJ);

      /* store 4 x x^T in lower half of JTJ */
      gsl_blas_dsyr(CblasLower, 4.0, x, JTJ);

      /* add alpha to diag(JTJ) */
      gsl_vector_add_constant(&diag.vector, par->alpha);
    }

  return GSL_SUCCESS;
}

int
penalty_fvv (const gsl_vector * x, const gsl_vector * v,
             void *params, gsl_vector * fvv)
{
  const size_t p = x->size;
  double normv = gsl_blas_dnrm2(v);

  gsl_vector_set_zero(fvv);
  gsl_vector_set(fvv, p, 2.0 * normv * normv);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

void
solve_system(const gsl_vector *x0, gsl_multilarge_nlinear_fdf *fdf,
             gsl_multilarge_nlinear_parameters *params)
{
  const gsl_multilarge_nlinear_type *T = gsl_multilarge_nlinear_trust;
  const size_t max_iter = 200;
  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 1.0e-8;
  const size_t n = fdf->n;
  const size_t p = fdf->p;
  gsl_multilarge_nlinear_workspace *work =
    gsl_multilarge_nlinear_alloc(T, params, n, p);
  gsl_vector * f = gsl_multilarge_nlinear_residual(work);
  gsl_vector * x = gsl_multilarge_nlinear_position(work);
  int info;
  double chisq0, chisq, rcond, xsq;
  struct timeval tv0, tv1;

  gettimeofday(&tv0, NULL);

  /* initialize solver */
  gsl_multilarge_nlinear_init(x0, fdf, work);

  /* store initial cost */
  gsl_blas_ddot(f, f, &chisq0);

  /* iterate until convergence */
  gsl_multilarge_nlinear_driver(max_iter, xtol, gtol, ftol,
                                NULL, NULL, &info, work);

  gettimeofday(&tv1, NULL);

  /* store final cost */
  gsl_blas_ddot(f, f, &chisq);

  /* compute final ||x||^2 */
  gsl_blas_ddot(x, x, &xsq);

  /* store cond(J(x)) */
  gsl_multilarge_nlinear_rcond(&rcond, work);

  /* print summary */
  fprintf(stderr, "%-25s %-5zu %-4zu %-5zu %-6zu %-4zu %-10.4e %-10.4e %-7.2f %-11.4e %.2f\n",
          gsl_multilarge_nlinear_trs_name(work),
          gsl_multilarge_nlinear_niter(work),
          fdf->nevalf,
          fdf->nevaldfu,
          fdf->nevaldf2,
          fdf->nevalfvv,
          chisq0,
          chisq,
          1.0 / rcond,
          xsq,
          (tv1.tv_sec - tv0.tv_sec) + 1.0e-6 * (tv1.tv_usec - tv0.tv_usec));

  gsl_multilarge_nlinear_free(work);
}

int
main (void)
{
  const size_t p = 2000;
  const size_t n = p + 1;
  gsl_vector *f = gsl_vector_alloc(n);
  gsl_vector *x = gsl_vector_alloc(p);

  /* allocate sparse Jacobian matrix with 2*p non-zero elements in triplet format */
  gsl_spmatrix *J = gsl_spmatrix_alloc_nzmax(n, p, 2 * p, GSL_SPMATRIX_TRIPLET);

  gsl_multilarge_nlinear_fdf fdf;
  gsl_multilarge_nlinear_parameters fdf_params =
    gsl_multilarge_nlinear_default_parameters();
  struct model_params params;
  size_t i;

  params.alpha = 1.0e-5;
  params.J = J;

  /* define function to be minimized */
  fdf.f = penalty_f;
  fdf.df = penalty_df;
  fdf.fvv = penalty_fvv;
  fdf.n = n;
  fdf.p = p;
  fdf.params = &params;

  for (i = 0; i < p; ++i)
    {
      /* starting point */
      gsl_vector_set(x, i, i + 1.0);

      /* store sqrt(alpha)*I_p in upper p-by-p block of J */
      gsl_spmatrix_set(J, i, i, sqrt(params.alpha));
    }

  fprintf(stderr, "%-25s %-4s %-4s %-5s %-6s %-4s %-10s %-10s %-7s %-11s %-10s\n",
          "Method", "NITER", "NFEV", "NJUEV", "NJTJEV", "NAEV", "Init Cost",
          "Final cost", "cond(J)", "Final |x|^2", "Time (s)");
  
  fdf_params.scale = gsl_multilarge_nlinear_scale_levenberg;

  fdf_params.trs = gsl_multilarge_nlinear_trs_lm;
  solve_system(x, &fdf, &fdf_params);

  fdf_params.trs = gsl_multilarge_nlinear_trs_lmaccel;
  solve_system(x, &fdf, &fdf_params);

  fdf_params.trs = gsl_multilarge_nlinear_trs_dogleg;
  solve_system(x, &fdf, &fdf_params);

  fdf_params.trs = gsl_multilarge_nlinear_trs_ddogleg;
  solve_system(x, &fdf, &fdf_params);

  fdf_params.trs = gsl_multilarge_nlinear_trs_subspace2D;
  solve_system(x, &fdf, &fdf_params);

  fdf_params.trs = gsl_multilarge_nlinear_trs_cgst;
  solve_system(x, &fdf, &fdf_params);

  gsl_vector_free(f);
  gsl_vector_free(x);
  gsl_spmatrix_free(J);

  return 0;
}
