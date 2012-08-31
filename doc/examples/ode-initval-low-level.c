int
main (void)
{
  const gsl_odeiv2_step_type * T 
    = gsl_odeiv2_step_rk8pd;

  gsl_odeiv2_step * s 
    = gsl_odeiv2_step_alloc (T, 2);
  gsl_odeiv2_control * c 
    = gsl_odeiv2_control_y_new (1e-6, 0.0);
  gsl_odeiv2_evolve * e 
    = gsl_odeiv2_evolve_alloc (2);

  double mu = 10;
  gsl_odeiv2_system sys = {func, jac, 2, &mu};

  double t = 0.0, t1 = 100.0;
  double h = 1e-6;
  double y[2] = { 1.0, 0.0 };

  while (t < t1)
    {
      int status = gsl_odeiv2_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, t1,
                                           &h, y);

      if (status != GSL_SUCCESS)
          break;

      printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  return 0;
}
