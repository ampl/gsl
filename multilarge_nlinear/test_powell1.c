#define powell1_N        4
#define powell1_P        4

static double powell1_x0[powell1_P] = { 3.0, -1.0, 0.0, 1.0 };
static double powell1_epsrel = 1.0e-4;

static double powell1_J[powell1_N * powell1_P];

static void
powell1_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < powell1_P; ++i)
    {
      gsl_test_rel(x[i], 0.0, epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
powell1_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get (x, 0);
  double x2 = gsl_vector_get (x, 1);
  double x3 = gsl_vector_get (x, 2);
  double x4 = gsl_vector_get (x, 3);

  gsl_vector_set(f, 0, x1 + 10.0*x2);
  gsl_vector_set(f, 1, sqrt(5.0) * (x3 - x4));
  gsl_vector_set(f, 2, pow(x2 - 2.0*x3, 2.0));
  gsl_vector_set(f, 3, sqrt(10.0) * pow((x1 - x4), 2.0));

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
powell1_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
            const gsl_vector * u, void * params, gsl_vector * v,
            gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(powell1_J, powell1_N, powell1_P);
  double x1 = gsl_vector_get (x, 0);
  double x2 = gsl_vector_get (x, 1);
  double x3 = gsl_vector_get (x, 2);
  double x4 = gsl_vector_get (x, 3);
  double term1 = x2 - 2.0*x3;
  double term2 = x1 - x4;

  gsl_matrix_set(&J.matrix, 0, 0, 1.0);
  gsl_matrix_set(&J.matrix, 0, 1, 10.0);
  gsl_matrix_set(&J.matrix, 0, 2, 0.0);
  gsl_matrix_set(&J.matrix, 0, 3, 0.0);

  gsl_matrix_set(&J.matrix, 1, 0, 0.0);
  gsl_matrix_set(&J.matrix, 1, 1, 0.0);
  gsl_matrix_set(&J.matrix, 1, 2, sqrt(5.0));
  gsl_matrix_set(&J.matrix, 1, 3, -sqrt(5.0));

  gsl_matrix_set(&J.matrix, 2, 0, 0.0);
  gsl_matrix_set(&J.matrix, 2, 1, 2.0*term1);
  gsl_matrix_set(&J.matrix, 2, 2, -4.0*term1);
  gsl_matrix_set(&J.matrix, 2, 3, 0.0);

  gsl_matrix_set(&J.matrix, 3, 0, 2.0*sqrt(10.0)*term2);
  gsl_matrix_set(&J.matrix, 3, 1, 0.0);
  gsl_matrix_set(&J.matrix, 3, 2, 0.0);
  gsl_matrix_set(&J.matrix, 3, 3, -2.0*sqrt(10.0)*term2);

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
powell1_fvv (const gsl_vector * x, const gsl_vector * v,
             void *params, gsl_vector * fvv)
{
  double v1 = gsl_vector_get(v, 0);
  double v2 = gsl_vector_get(v, 1);
  double v3 = gsl_vector_get(v, 2);
  double v4 = gsl_vector_get(v, 3);

  gsl_vector_set(fvv, 0, 0.0);
  gsl_vector_set(fvv, 1, 0.0);
  gsl_vector_set(fvv, 2, 2.0 * pow(v2 - 2.0*v3, 2.0));
  gsl_vector_set(fvv, 3, 2.0 * sqrt(10.0) * pow(v1 - v4, 2.0));

  (void)x;      /* avoid unused parameter warning */
  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf powell1_func =
{
  powell1_f,
  powell1_df,
  powell1_fvv,
  powell1_N,
  powell1_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem powell1_problem =
{
  "powell_singular",
  powell1_x0,
  NULL,
  &powell1_epsrel,
  &powell1_checksol,
  &powell1_func
};
