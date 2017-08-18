#define osborne_N         33
#define osborne_P         5

static double osborne_x0[osborne_P] = { 0.5, 1.5, -1.0, 0.01, 0.02 };
static double osborne_epsrel = 1.0e-8;

static double osborne_J[osborne_N * osborne_P];

static double osborne_Y[osborne_N] = {
0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881,
0.850, 0.818, 0.784, 0.751, 0.718, 0.685, 0.658,
0.628, 0.603, 0.580, 0.558, 0.538, 0.522, 0.506,
0.490, 0.478, 0.467, 0.457, 0.448, 0.438, 0.431,
0.424, 0.420, 0.414, 0.411, 0.406
};

static void
osborne_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  const double sumsq_exact = 5.464894697482687e-05;
  double osborne_x[osborne_P];

  osborne_x[0] = 3.754100521058740e-01;
  osborne_x[1] = GSL_NAN;
  osborne_x[2] = GSL_NAN;
  osborne_x[3] = GSL_NAN;
  osborne_x[4] = GSL_NAN;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  /* only the first model parameter is uniquely constrained */
  gsl_test_rel(x[0], osborne_x[0], epsrel, "%s/%s i=0",
               sname, pname);
}

static int
osborne_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  size_t i;

  for (i = 0; i < osborne_N; ++i)
    {
      double ti = 10.0*i;
      double yi = osborne_Y[i];
      double fi = yi - (x1 + x2*exp(-x4*ti) + x3*exp(-x5*ti));
      gsl_vector_set(f, i, fi);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
osborne_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
            const gsl_vector * u, void * params, gsl_vector * v,
            gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(osborne_J, osborne_N, osborne_P);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  size_t i;

  for (i = 0; i < osborne_N; ++i)
    {
      double ti = 10.0*i;
      double term1 = exp(-x4*ti);
      double term2 = exp(-x5*ti);

      gsl_matrix_set(&J.matrix, i, 0, -1.0);
      gsl_matrix_set(&J.matrix, i, 1, -term1);
      gsl_matrix_set(&J.matrix, i, 2, -term2);
      gsl_matrix_set(&J.matrix, i, 3, ti*x2*term1);
      gsl_matrix_set(&J.matrix, i, 4, ti*x3*term2);
    }

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
osborne_fvv (const gsl_vector * x, const gsl_vector * v,
             void *params, gsl_vector * fvv)
{
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double x4 = gsl_vector_get(x, 3);
  double x5 = gsl_vector_get(x, 4);
  double v2 = gsl_vector_get(v, 1);
  double v3 = gsl_vector_get(v, 2);
  double v4 = gsl_vector_get(v, 3);
  double v5 = gsl_vector_get(v, 4);
  size_t i;

  for (i = 0; i < osborne_N; ++i)
    {
      double ti = 10.0*i;
      double term1 = exp(-x4*ti);
      double term2 = exp(-x5*ti);
      double term3 = -2*v2 + ti*v4*x2;
      double term4 = -2*v3 + ti*v5*x3;

      gsl_vector_set(fvv, i, -term1 * term2 * ti *
                             (v4 / term2 * term3 + v5 / term1 * term4));
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf osborne_func =
{
  osborne_f,
  osborne_df,
  osborne_fvv,
  osborne_N,
  osborne_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem osborne_problem =
{
  "osborne",
  osborne_x0,
  NULL,
  &osborne_epsrel,
  &osborne_checksol,
  &osborne_func
};
