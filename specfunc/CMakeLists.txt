add_library(gslspecfunc
  airy.c airy_der.c airy_zero.c atanint.c bessel.c bessel.h bessel_I0.c
  bessel_I1.c bessel_In.c bessel_Inu.c bessel_J0.c bessel_J1.c bessel_Jn.c
  bessel_Jnu.c bessel_K0.c bessel_K1.c bessel_Kn.c bessel_Knu.c bessel_Y0.c
  bessel_Y1.c bessel_Yn.c bessel_Ynu.c bessel_amp_phase.c bessel_amp_phase.h
  bessel_i.c bessel_j.c bessel_k.c bessel_olver.c bessel_temme.c bessel_y.c
  bessel_zero.c bessel_sequence.c beta.c beta_inc.c clausen.c coulomb.c
  coupling.c coulomb_bound.c dawson.c debye.c dilog.c elementary.c ellint.c
  elljac.c erfc.c exp.c expint.c expint3.c fermi_dirac.c gegenbauer.c gamma.c
  gamma_inc.c hyperg_0F1.c hyperg_2F0.c hyperg_1F1.c hyperg_2F1.c hyperg_U.c
  hyperg.c laguerre.c lambert.c legendre_H3d.c legendre_Qn.c legendre_con.c
  legendre_poly.c log.c mathieu_angfunc.c mathieu_charv.c mathieu_coeff.c
  mathieu_radfunc.c mathieu_workspace.c poch.c pow_int.c psi.c recurse.h
  result.c shint.c sinint.c synchrotron.c transport.c trig.c zeta.c)

add_gsl_test(gslspecfunc-test
  test_sf.c test_sf.h test_airy.c test_bessel.c test_coulomb.c
  test_dilog.c test_gamma.c test_hyperg.c test_legendre.c test_mathieu.c
  LIBS gslspecfunc gsleigen gsllinalg gslsort gslmatrix gslvector gslblas
  gslcblas gslblock gslcomplex gslieeeutils gslerr gsltest gslsys utils)
