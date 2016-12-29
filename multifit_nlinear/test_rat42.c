#define rat42_N       9
#define rat42_P       3

static double rat42_x0a[rat42_P] = { 100.0, 1.0, 0.1 };
static double rat42_x0b[rat42_P] = { 75.0, 2.5, 0.07 };
static double rat42_epsrel = 1.0e-7;

static double rat42_sigma[rat42_P] = {
  1.7340283401E+00, 8.8295217536E-02, 3.4465663377E-03
};

static double rat42_X[rat42_N] = { 9.0, 14.0, 21.0, 28.0, 42.0,
                                   57.0, 63.0, 70.0, 79.0 };

static double rat42_F[rat42_N] = { 8.930, 10.800, 18.590, 22.330,
                                   39.350, 56.110, 61.730, 64.620,
                                   67.080 };

static void
rat42_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 8.0565229338E+00;
  const double rat42_x[rat42_P] = { 7.2462237576E+01,
                                    2.6180768402E+00,
                                    6.7359200066E-02 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < rat42_P; ++i)
    {
      gsl_test_rel(x[i], rat42_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}


static int
rat42_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double b[rat42_P];
  size_t i;

  for (i = 0; i < rat42_P; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < rat42_N; i++)
    {
      double xi = rat42_X[i];
      double yi = b[0] / (1.0 + exp(b[1] - b[2]*xi));
      gsl_vector_set (f, i, yi - rat42_F[i]);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
rat42_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  double b[rat42_P];
  size_t i;

  for (i = 0; i < rat42_P; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < rat42_N; i++)
    {
      double xi = rat42_X[i];
      double term1 = exp(b[1] - b[2]*xi);
      double term2 = 1.0 + term1;

      gsl_matrix_set (df, i, 0, 1.0 / term2);
      gsl_matrix_set (df, i, 1, -b[0] * term1 / (term2 * term2));
      gsl_matrix_set (df, i, 2, b[0] * term1 * xi / (term2 * term2));
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
rat42_fvv (const gsl_vector * x, const gsl_vector * v,
           void *params, gsl_vector * fvv)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double v1 = gsl_vector_get(v, 0);
  double v2 = gsl_vector_get(v, 1);
  double v3 = gsl_vector_get(v, 2);
  size_t i;

  for (i = 0; i < rat42_N; i++)
    {
      double ti = rat42_X[i];
      double term1 = exp(x2);
      double term2 = exp(ti * x3);

      gsl_vector_set(fvv, i,
        -pow(term1 + term2, -3.0) * term1 * term2 *
         (v2 - ti*v3) * (term1*(2*v1 - v2*x1 + ti*v3*x1) +
                         term2*(2*v1 + x1*(v2 - ti*v3))));
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multifit_nlinear_fdf rat42_func =
{
  rat42_f,
  rat42_df,
  rat42_fvv,
  rat42_N,
  rat42_P,
  NULL,
  0,
  0,
  0
};

static test_fdf_problem rat42a_problem =
{
  "nist-rat42a",
  rat42_x0a,
  NULL,
  rat42_sigma,
  &rat42_epsrel,
  &rat42_checksol,
  &rat42_func
};

static test_fdf_problem rat42b_problem =
{
  "nist-rat42b",
  rat42_x0b,
  NULL,
  rat42_sigma,
  &rat42_epsrel,
  &rat42_checksol,
  &rat42_func
};
