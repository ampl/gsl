#define powell2_N        2
#define powell2_P        2

static double powell2_x0[powell2_P] = { 3.0, 1.0 };
static double powell2_epsrel = 1.0e-3;

static double powell2_J[powell2_N * powell2_P];

static void
powell2_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < powell2_P; ++i)
    {
      gsl_test_rel(x[i], 0.0, epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
powell2_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);

  gsl_vector_set(f, 0, x0);
  gsl_vector_set(f, 1, 10.0*x0/(x0 + 0.1) + 2.0*x1*x1);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
powell2_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
            const gsl_vector * u, void * params, gsl_vector * v,
            gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(powell2_J, powell2_N, powell2_P);
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double term = x0 + 0.1;

  gsl_matrix_set(&J.matrix, 0, 0, 1.0);
  gsl_matrix_set(&J.matrix, 0, 1, 0.0);
  gsl_matrix_set(&J.matrix, 1, 0, 1.0 / (term * term));
  gsl_matrix_set(&J.matrix, 1, 1, 4.0 * x1);

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
powell2_fvv (const gsl_vector * x, const gsl_vector * v,
             void *params, gsl_vector * fvv)
{
  double x0 = gsl_vector_get (x, 0);
  double v0 = gsl_vector_get (v, 0);
  double v1 = gsl_vector_get (v, 1);
  double term = x0 + 0.1;

  gsl_vector_set(fvv, 0, 0.0);
  gsl_vector_set(fvv, 1, 4*v1*v1 - (2*v0*v0)/pow(term, 3.0));

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf powell2_func =
{
  powell2_f,
  powell2_df,
  powell2_fvv,
  powell2_N,
  powell2_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem powell2_problem =
{
  "powell2",
  powell2_x0,
  NULL,
  &powell2_epsrel,
  &powell2_checksol,
  &powell2_func
};
