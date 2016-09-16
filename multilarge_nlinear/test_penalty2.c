#define penalty2_N         8 /* 2*p */
#define penalty2_P         4

static double penalty2_x0[penalty2_P] = { 0.5, 0.5, 0.5, 0.5 };
static double penalty2_epsrel = 1.0e-12;

static double penalty2_J[penalty2_N * penalty2_P];

static void
penalty2_checksol(const double x[], const double sumsq,
                  const double epsrel, const char *sname,
                  const char *pname)
{
  const double sumsq_exact = 9.37629300735544219e-06;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  (void)x; /* avoid unused parameter warning */
}

static int
penalty2_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  const double alpha = 1.0e-5;
  const double sqrt_alpha = sqrt(alpha);
  double x1 = gsl_vector_get(x, 0);
  size_t i;
  double sum = penalty2_P * x1 * x1;

  gsl_vector_set(f, 0, x1 - 0.2);

  /* rows [2:p] */
  for (i = 1; i < penalty2_P; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double xim1 = gsl_vector_get(x, i - 1);
      double yi = exp(0.1*(i + 1.0)) + exp(0.1*i);

      gsl_vector_set(f, i, sqrt_alpha*(exp(0.1*xi) + exp(0.1*xim1) - yi));

      sum += (penalty2_P - i) * xi * xi;
    }

  /* rows [p+1:2p-1] */
  for (i = penalty2_P; i < penalty2_N - 1; ++i)
    {
      double xi = gsl_vector_get(x, i - penalty2_P + 1);

      gsl_vector_set(f, i, sqrt_alpha*(exp(0.1*xi) - exp(-0.1)));
    }

  /* row 2p */
  gsl_vector_set(f, penalty2_N - 1, sum - 1.0);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
penalty2_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
             const gsl_vector * u, void * params, gsl_vector * v,
             gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(penalty2_J, penalty2_N, penalty2_P);
  const double alpha = 1.0e-5;
  const double sqrt_alpha = sqrt(alpha);
  size_t i, j;

  for (j = 0; j < penalty2_P; ++j)
    {
      double xj = gsl_vector_get(x, j);
      double delta1j = (j == 0) ? 1.0 : 0.0;

      /* first and last rows */
      gsl_matrix_set(&J.matrix, 0, j, delta1j);
      gsl_matrix_set(&J.matrix, penalty2_N - 1, j, 2.0 * (penalty2_P - j) * xj);

      /* rows [2:p] */
      for (i = 1; i < penalty2_P; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double xim1 = gsl_vector_get(x, i - 1);
          double Jij;

          if (i == j)
            Jij = exp(0.1 * xi);
          else if (i - 1 == j)
            Jij = exp(0.1 * xim1);
          else
            Jij = 0.0;

          Jij *= 0.1 * sqrt_alpha;

          gsl_matrix_set(&J.matrix, i, j, Jij);
        }

      /* rows [p+1:2p-1] */
      for (i = penalty2_P; i < penalty2_N - 1; ++i)
        {
          double xi = gsl_vector_get(x, i - penalty2_P + 1);

          if (i - penalty2_P + 1 == j)
            gsl_matrix_set(&J.matrix, i, j, 0.1 * sqrt_alpha * exp(0.1*xi));
          else
            gsl_matrix_set(&J.matrix, i, j, 0.0);
        }
    }

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
penalty2_fvv (const gsl_vector * x, const gsl_vector * v,
              void *params, gsl_vector * fvv)
{
  const double alpha = 1.0e-5;
  const double sqrt_alpha = sqrt(alpha);
  double v1 = gsl_vector_get(v, 0);
  double sum = penalty2_P * v1 * v1;
  size_t i;

  /* first row */
  gsl_vector_set(fvv, 0, 0.0);

  /* rows [2:p] */
  for (i = 1; i < penalty2_P; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double xim1 = gsl_vector_get(x, i - 1);
      double vi = gsl_vector_get(v, i);
      double vim1 = gsl_vector_get(v, i - 1);
      double term1 = exp(xi / 10.0);
      double term2 = exp(xim1 / 10.0);

      gsl_vector_set(fvv, i,
        sqrt_alpha / 100.0 * (term1 * vi * vi + term2 * vim1 * vim1));

      sum += (penalty2_P - i) * vi * vi;
    }

  /* last row */
  gsl_vector_set(fvv, penalty2_N - 1, 2.0 * sum);

  /* rows [p+1:2p-1] */
  for (i = penalty2_P; i < penalty2_N - 1; ++i)
    {
      double xi = gsl_vector_get(x, i - penalty2_P + 1);
      double vi = gsl_vector_get(v, i - penalty2_P + 1);

      gsl_vector_set(fvv, i, sqrt_alpha / 100.0 * exp(xi / 10.0) * vi * vi);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf penalty2_func =
{
  penalty2_f,
  penalty2_df,
  penalty2_fvv,
  penalty2_N,
  penalty2_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem penalty2_problem =
{
  "penalty2",
  penalty2_x0,
  NULL,
  &penalty2_epsrel,
  &penalty2_checksol,
  &penalty2_func
};
