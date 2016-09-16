#define watson_N         31
#define watson_P         6

static double watson_x0[watson_P] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double watson_epsrel = 1.0e-6;

static double watson_J[watson_N * watson_P];

static void
watson_checksol(const double x[], const double sumsq,
                const double epsrel, const char *sname,
                const char *pname)
{
  size_t i;
  const double sumsq_exact = 2.287670053552372e-03;
  const double watson_x[watson_P] = {
    -1.572508640629858e-02,  1.012434869366059e+00, -2.329916259263380e-01,
     1.260430087686035e+00, -1.513728922580576e+00,  9.929964323646112e-01 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < watson_P; ++i)
    {
      gsl_test_rel(x[i], watson_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
watson_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  size_t i, j;

  for (i = 0; i < watson_N - 2; ++i)
    {
      double ti = (i + 1) / 29.0;
      double tjm1 = 1.0, tjm2 = 1.0;
      double sum1 = 0.0, sum2 = 0.0;

      for (j = 0; j < watson_P; ++j)
        {
          double xj = gsl_vector_get(x, j);

          sum1 += xj * tjm1;
          tjm1 *= ti;

          if (j > 0)
            {
              sum2 += j * xj * tjm2;
              tjm2 *= ti;
            }
        }

      gsl_vector_set (f, i, sum2 - sum1*sum1 - 1.0);
    }

  gsl_vector_set(f, watson_N - 2, x1);
  gsl_vector_set(f, watson_N - 1, x2 - x1*x1 - 1.0);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
watson_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
           const gsl_vector * u, void * params, gsl_vector * v,
           gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(watson_J, watson_N, watson_P);
  double x1 = gsl_vector_get (x, 0);
  size_t i, j;

  gsl_matrix_set_zero(&J.matrix);

  for (i = 0; i < watson_N - 2; ++i)
    {
      double ti = (i + 1) / 29.0;
      double tjm1 = 1.0, tjm2 = 1.0;
      double sum1 = 0.0;

      for (j = 0; j < watson_P; ++j)
        {
          double xj = gsl_vector_get(x, j);
          sum1 += xj * tjm1;
          tjm1 *= ti;
        }

      tjm1 = 1.0;
      tjm2 = 1.0;
      for (j = 0; j < watson_P; ++j)
        {
          gsl_matrix_set(&J.matrix, i, j, j * tjm2 - 2.0*sum1*tjm1);
          tjm1 *= ti;

          if (j > 0)
            tjm2 *= ti;
        }
    }

  gsl_matrix_set(&J.matrix, watson_N - 2, 0, 1.0);
  gsl_matrix_set(&J.matrix, watson_N - 1, 0, -2.0*x1);
  gsl_matrix_set(&J.matrix, watson_N - 1, 1, 1.0);

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
watson_fvv (const gsl_vector * x, const gsl_vector * v,
            void *params, gsl_vector * fvv)
{
  double v1 = gsl_vector_get (v, 0);
  size_t i, j;

  for (i = 0; i < watson_N - 2; ++i)
    {
      double ti = (i + 1) / 29.0;
      double sum = 0.0;
      double tjm1 = 1.0;

      for (j = 0; j < watson_P; ++j)
        {
          double vj = gsl_vector_get(v, j);
          sum += vj * tjm1;
          tjm1 *= ti;
        }

      gsl_vector_set(fvv, i, -2.0*sum*sum);
    }

  gsl_vector_set(fvv, watson_N - 2, 0.0);
  gsl_vector_set(fvv, watson_N - 1, -2.0*v1*v1);

  (void)x;      /* avoid unused parameter warning */
  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf watson_func =
{
  watson_f,
  watson_df,
  watson_fvv,
  watson_N,
  watson_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem watson_problem =
{
  "watson",
  watson_x0,
  NULL,
  &watson_epsrel,
  &watson_checksol,
  &watson_func
};
