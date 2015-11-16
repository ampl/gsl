/* expfit.c -- model functions for exponential + background */

struct data {
  size_t n;
  double * y;
};

int
expb_f (const gsl_vector * x, void *data, 
        gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double b = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      double t = i;
      double Yi = A * exp (-lambda * t) + b;
      gsl_vector_set (f, i, Yi - y[i]);
    }

  return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *data, 
         gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      double t = i;
      double e = exp(-lambda * t);
      gsl_matrix_set (J, i, 0, e); 
      gsl_matrix_set (J, i, 1, -t * A * e);
      gsl_matrix_set (J, i, 2, 1.0);
    }
  return GSL_SUCCESS;
}
