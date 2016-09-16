#define boxbod_N       6
#define boxbod_P       2

static double boxbod_x0a[boxbod_P] = { 1.0, 1.0 };
static double boxbod_x0b[boxbod_P] = { 100.0, 0.75 };
static double boxbod_epsrel = 1.0e-7;

static double boxbod_J[boxbod_N * boxbod_P];

static double boxbod_sigma[boxbod_P] = {
  1.2354515176E+01, 1.0455993237E-01
};

static double boxbod_X[boxbod_N] = { 1.0, 2.0, 3.0, 5.0, 7.0, 10.0 };

static double boxbod_F[boxbod_N] = { 109.0, 149.0, 149.0, 191.0,
                                     213.0, 224.0 };

static void
boxbod_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 1.1680088766E+03;
  const double boxbod_x[boxbod_P] = { 2.1380940889E+02,
                                      5.4723748542E-01 };

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < boxbod_P; ++i)
    {
      gsl_test_rel(x[i], boxbod_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}


static int
boxbod_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double b[boxbod_P];
  size_t i;

  for (i = 0; i < boxbod_P; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < boxbod_N; i++)
    {
      double xi = boxbod_X[i];
      double yi;

      yi = b[0] * (1.0 - exp(-b[1] * xi));
      gsl_vector_set (f, i, yi - boxbod_F[i]);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
boxbod_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
           const gsl_vector * u, void * params, gsl_vector * v,
           gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(boxbod_J, boxbod_N, boxbod_P);
  double b[boxbod_P];
  size_t i;

  for (i = 0; i < boxbod_P; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < boxbod_N; i++)
    {
      double xi = boxbod_X[i];
      double term = exp(-b[1] * xi);

      gsl_matrix_set (&J.matrix, i, 0, 1.0 - term);
      gsl_matrix_set (&J.matrix, i, 1, b[0] * term * xi);
    }

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
boxbod_fvv (const gsl_vector * x, const gsl_vector * v,
            void *params, gsl_vector * fvv)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double v1 = gsl_vector_get(v, 0);
  double v2 = gsl_vector_get(v, 1);
  size_t i;

  for (i = 0; i < boxbod_N; i++)
    {
      double ti = boxbod_X[i];
      double term = exp(-x2 * ti);

      gsl_vector_set(fvv, i, term * ti * v2 * (2*v1 - ti*v2*x1));
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf boxbod_func =
{
  boxbod_f,
  boxbod_df,
  boxbod_fvv,
  boxbod_N,
  boxbod_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem boxboda_problem =
{
  "nist-boxboda",
  boxbod_x0a,
  boxbod_sigma,
  &boxbod_epsrel,
  &boxbod_checksol,
  &boxbod_func
};

static test_fdf_problem boxbodb_problem =
{
  "nist-boxbodb",
  boxbod_x0b,
  boxbod_sigma,
  &boxbod_epsrel,
  &boxbod_checksol,
  &boxbod_func
};
