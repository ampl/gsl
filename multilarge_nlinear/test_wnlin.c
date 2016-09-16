#define wnlin_N         40
#define wnlin_P         3

static double wnlin_x0[wnlin_P] = { 1.0, 0.9, 0.0 };
static double wnlin_epsrel = 1.0e-8;

static double wnlin_J[wnlin_N * wnlin_P];

/* data */
static double wnlin_Y[wnlin_N] = {
  6.08035e+00, 5.47552e+00, 5.94654e+00, 5.04920e+00, 4.78568e+00,
  3.51748e+00, 2.84671e+00, 3.24634e+00, 3.23395e+00, 3.30385e+00,
  2.83439e+00, 2.31891e+00, 2.33858e+00, 2.40559e+00, 2.41856e+00,
  1.99966e+00, 1.88127e+00, 1.91477e+00, 1.70415e+00, 1.60316e+00,
  1.77937e+00, 1.55302e+00, 1.50903e+00, 1.36364e+00, 1.36873e+00,
  1.41954e+00, 1.37778e+00, 1.23573e+00, 1.28524e+00, 1.46327e+00,
  1.22315e+00, 1.19330e+00, 1.18717e+00, 8.83172e-01, 1.23424e+00,
  1.14683e+00, 1.11091e+00, 1.20396e+00, 1.28722e+00, 1.05801e+00
};

/* weights */
static double wnlin_W[wnlin_N] = {
  2.77778e+00, 3.27690e+00, 3.85426e+00, 4.51906e+00, 5.28083e+00,
  6.14919e+00, 7.13370e+00, 8.24349e+00, 9.48703e+00, 1.08717e+01,
  1.24036e+01, 1.40869e+01, 1.59238e+01, 1.79142e+01, 2.00553e+01,
  2.23415e+01, 2.47646e+01, 2.73137e+01, 2.99753e+01, 3.27337e+01,
  3.55714e+01, 3.84696e+01, 4.14085e+01, 4.43678e+01, 4.73278e+01,
  5.02690e+01, 5.31731e+01, 5.60234e+01, 5.88046e+01, 6.15036e+01,
  6.41092e+01, 6.66121e+01, 6.90054e+01, 7.12839e+01, 7.34442e+01,
  7.54848e+01, 7.74053e+01, 7.92069e+01, 8.08918e+01, 8.24632e+01
};

static void
wnlin_checksol(const double x[], const double sumsq,
               const double epsrel, const char *sname,
               const char *pname)
{
  size_t i;
  const double sumsq_exact = 29.7481259665713758;
  const double wnlin_x[wnlin_P] = { 5.17378551196259195,
                                    0.111041758006851149,
                                    1.05282724070446099 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < wnlin_P; ++i)
    {
      gsl_test_rel(x[i], wnlin_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
wnlin_f (const gsl_vector *x, void *params, gsl_vector *f)
{
  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double b = gsl_vector_get (x, 2);
  size_t i;

  /* model Yi = A * exp(-lambda * i) + b */
  for (i = 0; i < wnlin_N; i++)
    {
      double ti = i;
      double yi = wnlin_Y[i];
      double swi = sqrt(wnlin_W[i]);
      double Mi = A * exp (-lambda * ti) + b;

      gsl_vector_set (f, i, swi * (Mi - yi));
    }

  return GSL_SUCCESS;
}

static int
wnlin_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
          const gsl_vector * u, void * params, gsl_vector * v,
          gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(wnlin_J, wnlin_N, wnlin_P);
  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  size_t i;

  for (i = 0; i < wnlin_N; i++)
    {
      gsl_vector_view v = gsl_matrix_row(&J.matrix, i);
      double ti = i;
      double swi = sqrt(wnlin_W[i]);
      double e = exp(-lambda * ti);

      gsl_vector_set(&v.vector, 0, e);
      gsl_vector_set(&v.vector, 1, -ti * A * e);
      gsl_vector_set(&v.vector, 2, 1.0);

      gsl_vector_scale(&v.vector, swi);
    }

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  return GSL_SUCCESS;
}

static int
wnlin_fvv (const gsl_vector * x, const gsl_vector * v,
           void *params, gsl_vector * fvv)
{
  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double v1 = gsl_vector_get(v, 0);
  double v2 = gsl_vector_get(v, 1);
  size_t i;

  for (i = 0; i < wnlin_N; i++)
    {
      double ti = i;
      double swi = sqrt(wnlin_W[i]);
      double fvvi;

      fvvi = exp(-ti*lambda)*ti*v2 * (-2*v1 + ti*v2*A);

      gsl_vector_set(fvv, i, swi * fvvi);
    }

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf wnlin_func1 =
{
  wnlin_f,
  wnlin_df,
  wnlin_fvv,
  wnlin_N,
  wnlin_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem wnlin_problem1 =
{
  "wnlin_internal_weights",
  wnlin_x0,
  NULL,
  &wnlin_epsrel,
  &wnlin_checksol,
  &wnlin_func1
};
