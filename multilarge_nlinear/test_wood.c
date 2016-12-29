#define wood_N         6
#define wood_P         4

static double wood_x0[wood_P] = { -3.0, -1.0, -3.0, -1.0 };
static double wood_epsrel = 1.0e-12;

static double wood_J[wood_N * wood_P];

static void
wood_checksol(const double x[], const double sumsq,
              const double epsrel, const char *sname,
              const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;
  const double wood_x[wood_P] = { 1.0, 1.0, 1.0, 1.0 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < wood_P; ++i)
    {
      gsl_test_rel(x[i], wood_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
wood_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);

  gsl_vector_set(f, 0, 10.0*(x2 - x1*x1));
  gsl_vector_set(f, 1, 1.0 - x1);
  gsl_vector_set(f, 2, sqrt(90.0)*(x4 - x3*x3));
  gsl_vector_set(f, 3, 1.0 - x3);
  gsl_vector_set(f, 4, sqrt(10.0)*(x2 + x4 - 2.0));
  gsl_vector_set(f, 5, (x2 - x4) / sqrt(10.0));

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
wood_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
         const gsl_vector * u, void * params, gsl_vector * v,
         gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(wood_J, wood_N, wood_P);
  double x1 = gsl_vector_get(x, 0);
  double x3 = gsl_vector_get(x, 2);
  double s90 = sqrt(90.0);
  double s10 = sqrt(10.0);

  gsl_matrix_set_zero(&J.matrix);

  gsl_matrix_set(&J.matrix, 0, 0, -20.0*x1);
  gsl_matrix_set(&J.matrix, 0, 1, 10.0);
  gsl_matrix_set(&J.matrix, 1, 0, -1.0);
  gsl_matrix_set(&J.matrix, 2, 2, -2.0*s90*x3);
  gsl_matrix_set(&J.matrix, 2, 3, s90);
  gsl_matrix_set(&J.matrix, 3, 2, -1.0);
  gsl_matrix_set(&J.matrix, 4, 1, s10);
  gsl_matrix_set(&J.matrix, 4, 3, s10);
  gsl_matrix_set(&J.matrix, 5, 1, 1.0/s10);
  gsl_matrix_set(&J.matrix, 5, 3, -1.0/s10);

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
wood_fvv (const gsl_vector * x, const gsl_vector * v,
          void *params, gsl_vector * fvv)
{
  const double s10 = sqrt(10.0);
  double v1 = gsl_vector_get(v, 0);
  double v3 = gsl_vector_get(v, 2);

  gsl_vector_set(fvv, 0, -20.0 * v1 * v1);
  gsl_vector_set(fvv, 1, 0.0);
  gsl_vector_set(fvv, 2, -6.0 * s10 * v3 * v3);
  gsl_vector_set(fvv, 3, 0.0);
  gsl_vector_set(fvv, 4, 0.0);
  gsl_vector_set(fvv, 5, 0.0);

  (void)x;      /* avoid unused parameter warning */
  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf wood_func =
{
  wood_f,
  wood_df,
  wood_fvv,
  wood_N,
  wood_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem wood_problem =
{
  "wood",
  wood_x0,
  NULL,
  &wood_epsrel,
  &wood_checksol,
  &wood_func
};
