int
main (void)
{
  double mu = 10;
  gsl_odeiv2_system sys = { func, jac, 2, &mu };

  gsl_odeiv2_driver *d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
                                   1e-3, 1e-8, 1e-8);

  double t = 0.0;
  double y[2] = { 1.0, 0.0 };
  int i, s;

  for (i = 0; i < 100; i++)
    {
      s = gsl_odeiv2_driver_apply_fixed_step (d, &t, 1e-3, 1000, y);

      if (s != GSL_SUCCESS)
        {
          printf ("error: driver returned %d\n", s);
          break;
        }

      printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

  gsl_odeiv2_driver_free (d);
  return s;
}
