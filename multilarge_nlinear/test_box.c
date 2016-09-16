#define box_N         10 /* can be >= p */
#define box_P         3

/* dogleg method fails with recommended starting point, so use
 * a slightly easier x0 */
/*static double box_x0[box_P] = { 0.0, 10.0, 20.0 };*/
static double box_x0[box_P] = { 5.0, 10.0, 2.0 };
static double box_epsrel = 1.0e-12;

static double box_J[box_N * box_P];

static void
box_checksol(const double x[], const double sumsq,
              const double epsrel, const char *sname,
              const char *pname)
{
  const double sumsq_exact = 0.0;
  const double eps = 1.0e-6;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  /* there are 3 possible solution vectors */

  if (fabs(x[2] - 1.0) < eps)
    {
      /* case 1: x* = [ 1; 10; 1 ] */
      gsl_test_rel(x[0], 1.0, epsrel, "%s/%s i=0",
                   sname, pname);
      gsl_test_rel(x[1], 10.0, epsrel, "%s/%s i=1",
                   sname, pname);
      gsl_test_rel(x[2], 1.0, epsrel, "%s/%s i=2",
                   sname, pname);
    }
  else if (fabs(x[2] + 1.0) < eps)
    {
      /* case 2: x* = [ 10; 1; -1 ] */
      gsl_test_rel(x[0], 10.0, epsrel, "%s/%s i=0",
                   sname, pname);
      gsl_test_rel(x[1], 1.0, epsrel, "%s/%s i=1",
                   sname, pname);
      gsl_test_rel(x[2], -1.0, epsrel, "%s/%s i=2",
                   sname, pname);
    }
  else
    {
      /* case 3: x* = [ a; a; 0 ] for any a */
      gsl_test_rel(x[0], x[1], epsrel, "%s/%s i=0,1",
                   sname, pname);
      gsl_test_rel(x[2], 0.0, epsrel, "%s/%s i=2",
                   sname, pname);
    }
}

static int
box_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  size_t i;

  for (i = 0; i < box_N; ++i)
    {
      double ti = (i + 1.0) / 10.0;
      double fi = exp(-x1*ti) - exp(-x2*ti) - x3*(exp(-ti) - exp(-10.0*ti));
      gsl_vector_set(f, i, fi);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
box_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
        const gsl_vector * u, void * params, gsl_vector * v,
        gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(box_J, box_N, box_P);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  size_t i;

  for (i = 0; i < box_N; ++i)
    {
      double ti = (i + 1.0) / 10.0;
      double term1 = exp(-x1*ti);
      double term2 = exp(-x2*ti);
      double term3 = exp(-10.0*ti) - exp(-ti);

      gsl_matrix_set(&J.matrix, i, 0, -ti*term1);
      gsl_matrix_set(&J.matrix, i, 1, ti*term2);
      gsl_matrix_set(&J.matrix, i, 2, term3);
    }

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
box_fvv (const gsl_vector * x, const gsl_vector * v,
         void *params, gsl_vector * fvv)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double v1 = gsl_vector_get(v, 0);
  double v2 = gsl_vector_get(v, 1);
  size_t i;

  for (i = 0; i < box_N; ++i)
    {
      double ti = (i + 1.0) / 10.0;
      double term1 = exp(-x1*ti);
      double term2 = exp(-x2*ti);

      gsl_vector_set(fvv, i, ti * ti * (v1*v1*term1 - v2*v2*term2));
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf box_func =
{
  box_f,
  box_df,
  box_fvv,
  box_N,
  box_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem box_problem =
{
  "box3d",
  box_x0,
  NULL,
  &box_epsrel,
  &box_checksol,
  &box_func
};
