#define biggs_N         20 /* >= p */
#define biggs_P         6

/* dogleg method has trouble converging from recommended starting point,
 * so we use an x0 which is a little closer to the true solution */
/*static double biggs_x0[biggs_P] = { 1.0, 2.0, 1.0, 1.0, 1.0, 1.0 };*/
static double biggs_x0[biggs_P] = { 1.0, 5.0, 1.0, 2.0, 3.0, 2.0 };
static double biggs_epsrel = 1.0e-9;

static void
biggs_checksol(const double x[], const double sumsq,
               const double epsrel, const char *sname,
               const char *pname)
{
  const double sumsq_exact = 0.0;
  const double biggs_x[biggs_P] = { 1.0, 10.0, 1.0, 5.0, 4.0, 3.0 };
  size_t i;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < biggs_P; ++i)
    {
      gsl_test_rel(x[i], biggs_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
biggs_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  double x6 = gsl_vector_get(x, 5);
  size_t i;

  for (i = 0; i < biggs_N; ++i)
    {
      double ti = 0.1 * (i + 1.0);
      double yi = exp(-ti) - 5*exp(-10*ti) + 3*exp(-4*ti);
      double fi = x3*exp(-ti*x1) - x4*exp(-ti*x2) + x6*exp(-ti*x5) - yi;

      gsl_vector_set(f, i, fi);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
biggs_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  double x6 = gsl_vector_get(x, 5);
  size_t i;

  for (i = 0; i < biggs_N; ++i)
    {
      double ti = 0.1 * (i + 1.0);

      gsl_matrix_set(J, i, 0, -ti*x3*exp(-ti*x1));
      gsl_matrix_set(J, i, 1, ti*x4*exp(-ti*x2));
      gsl_matrix_set(J, i, 2, exp(-ti*x1));
      gsl_matrix_set(J, i, 3, -exp(-ti*x2));
      gsl_matrix_set(J, i, 4, -ti*x6*exp(-ti*x5));
      gsl_matrix_set(J, i, 5, exp(-ti*x5));
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
biggs_fvv (const gsl_vector * x, const gsl_vector * v,
           void *params, gsl_vector * fvv)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  double x6 = gsl_vector_get(x, 5);
  double v1 = gsl_vector_get(v, 0);
  double v2 = gsl_vector_get(v, 1);
  double v3 = gsl_vector_get(v, 2);
  double v4 = gsl_vector_get(v, 3);
  double v5 = gsl_vector_get(v, 4);
  double v6 = gsl_vector_get(v, 5);
  size_t i;

  for (i = 0; i < biggs_N; ++i)
    {
      double ti = 0.1 * (i + 1.0);
      double term1 = exp(-ti * x1);
      double term2 = exp(-ti * x2);
      double term3 = exp(-ti * x5);

      gsl_vector_set(fvv, i, ti * term1 * term2 * term3 *
                             (v1/(term2*term3)*(-2*v3 + ti*v1*x3) -
                              v2/(term1*term3)*(-2*v4 + ti*v2*x4) +
                              v5/(term1*term2)*(-2*v6 + ti*v5*x6)));
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multifit_nlinear_fdf biggs_func =
{
  biggs_f,
  biggs_df,
  biggs_fvv,
  biggs_N,
  biggs_P,
  NULL,
  0,
  0,
  0
};

static test_fdf_problem biggs_problem =
{
  "biggs",
  biggs_x0,
  NULL,
  NULL,
  &biggs_epsrel,
  &biggs_checksol,
  &biggs_func
};
