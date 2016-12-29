#define brown1_N         20
#define brown1_P         4

static double brown1_x0[brown1_P] = { 25, 5, -5, -1 };
static double brown1_epsrel = 1.0e-5;

static void
brown1_checksol(const double x[], const double sumsq,
                const double epsrel, const char *sname,
                const char *pname)
{
  size_t i;
  const double sumsq_exact = 8.582220162635628e+04;
  const double brown1_x[brown1_P] = {
    -1.159443990239263e+01, 1.320363005221244e+01,
    -4.034395456782477e-01, 2.367789088597534e-01 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < brown1_P; ++i)
    {
      gsl_test_rel(x[i], brown1_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
brown1_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double x2 = gsl_vector_get (x, 2);
  double x3 = gsl_vector_get (x, 3);
  size_t i;

  for (i = 0; i < brown1_N; i++)
    {
      double ti = 0.2 * (i + 1);
      double ui = x0 + x1 * ti - exp (ti);
      double vi = x2 + x3 * sin (ti) - cos (ti);

      gsl_vector_set (f, i, ui * ui + vi * vi);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
brown1_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double x2 = gsl_vector_get (x, 2);
  double x3 = gsl_vector_get (x, 3);
  size_t i;

  for (i = 0; i < brown1_N; i++)
    {
      double ti = 0.2 * (i + 1);
      double ui = x0 + x1 * ti - exp (ti);
      double vi = x2 + x3 * sin (ti) - cos (ti);

      gsl_matrix_set (df, i, 0, 2 * ui);
      gsl_matrix_set (df, i, 1, 2 * ui * ti);
      gsl_matrix_set (df, i, 2, 2 * vi);
      gsl_matrix_set (df, i, 3, 2 * vi * sin (ti));

    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
brown1_fvv (const gsl_vector * x, const gsl_vector * v,
            void *params, gsl_vector * fvv)
{
  double v0 = gsl_vector_get (v, 0);
  double v1 = gsl_vector_get (v, 1);
  double v2 = gsl_vector_get (v, 2);
  double v3 = gsl_vector_get (v, 3);
  size_t i;

  for (i = 0; i < brown1_N; i++)
    {
      double ti = 0.2 * (i + 1);
      double term1 = v0 + ti*v1;
      double term2 = v3*sin(ti);

      gsl_vector_set (fvv, i, 2.0 * (term1*term1 + v2*v2 +
                                     term2 * (2*v2 + term2)));
    }

  (void)x;      /* avoid unused parameter warning */
  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multifit_nlinear_fdf brown1_func =
{
  brown1_f,
  brown1_df,
  brown1_fvv,
  brown1_N,
  brown1_P,
  NULL,
  0,
  0,
  0
};

static test_fdf_problem brown1_problem =
{
  "brown_dennis",
  brown1_x0,
  NULL,
  NULL,
  &brown1_epsrel,
  &brown1_checksol,
  &brown1_func
};
