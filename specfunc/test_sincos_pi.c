/* specfunc/test_sincos_pi.c
 * 
 * Copyright (C) 2017 Konrad Griessinger
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author: Konrad Griessinger */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf.h>
#include "test_sf.h"

/* Any double precision number bigger than this is automatically an even integer. */
#define BIGDBL (2.0 / GSL_DBL_EPSILON)

int
test_sincos_pi(void)
{
  gsl_sf_result r;
  int s = 0;
  int k = 0, kmax = 12;
  double x = 0.0, ix = 0.0, fx = 0.0, exact = 0.0;

  /* sin_pi tests */

  fx = 0.5;
  exact = 1.0;

  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
  fx = -0.5;
  exact = -1.0;

  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = 1.5;
  exact = -1.0;

  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
  fx = -1.5;
  exact = 1.0;

  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = 2.5;
  exact = 1.0;

  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
  fx = -2.5;
  exact = -1.0;

  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = 3.5;
  exact = -1.0;

  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
  fx = -3.5;
  exact = 1.0;

  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = 0.375;
  exact = 0.923879532511286756128183189397;
  
  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = -0.375;
  exact = -0.923879532511286756128183189397;
  
  TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
  
  fx = 0.0;
  exact = 0.0;

  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  fx = 0.5;
  exact = 1.0;

  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  fx = 0.03125;
  exact = 0.0980171403295606019941955638886;

  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  
  fx = 0.0625;
  exact = 0.195090322016128267848284868477;

  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }

  fx = 0.75;
  exact = 0.707106781186547524400844362105;

  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  fx = 0.0078125;
  exact = 0.0245412285229122880317345294593;
  
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }

  /* sin_pi tests for very large arguments */

  fx = 0.0625;
  exact = 0.195090322016128267848284868477;
  ix = LONG_MAX + 1.0;
  ix += fabs(fmod(ix,2.0)); /* make sure of even number */
  
  for (k=0; k<kmax; k++) {
    x = ix + fx;
    x -= ix; /* careful with compiler optimization */
    if ( ( x != fx ) || ( fabs(ix+fx) >= BIGDBL ) ) break;
    printf("ix+fx= %.18e\n", ix+fx);
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix += 101.0;
    exact = -exact;
  }

  fx = -0.0625;
  exact = -0.195090322016128267848284868477;
  ix = LONG_MIN - 1.0;
  ix -= fabs(fmod(ix,2.0)); /* make sure of even number */
  
  for (k=0; k<kmax; k++) {
    x = ix + fx;
    x -= ix; /* careful with compiler optimization */
    if ( ( x != fx ) || ( fabs(ix+fx) >= BIGDBL ) ) break;
    printf("ix+fx= %.18e\n", ix+fx);
    TEST_SF(s, gsl_sf_sin_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix -= 101.0;
    exact = -exact;
  }

  
  
  /* cos_pi tests */

  ix = 0.0;
  fx = 0.0;
  exact = 1.0;

  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = 1.0;
  exact = -1.0;

  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = -1.0;
  exact = -1.0;

  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = 2.0;
  exact = 1.0;

  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = -2.0;
  exact = 1.0;

  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = 3.0;
  exact = -1.0;

  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = -3.0;
  exact = -1.0;

  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);


  fx = 0.375;
  exact = 0.382683432365089771728459984030;
  
  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);

  fx = -0.375;
  exact = 0.382683432365089771728459984030;
  
  TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
  
  fx = 0.0;
  exact = 1.0;
  
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  fx = 0.5;
  exact = 0.0;
  
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  fx = 0.0625;
  exact = 0.980785280403230449126182236134;
  
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  fx = 0.4375;
  exact = 0.195090322016128267848284868477;
  
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  fx = 0.4921875;
  exact = 0.0245412285229122880317345294593;
  
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(3.0,k+1);
    if (k==0) exact = -exact;
  }

  exact = fabs(exact);
  ix = 0.0;
  for (k=0; k<kmax; k++) {
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix = pow(10.0,k+1);
  }
  
  /* cos_pi tests for very large arguments */

  fx = 0.0625;
  exact = 0.980785280403230449126182236134;
  ix = LONG_MAX + 1.0;
  ix += fabs(fmod(ix,2.0)); /* make sure of even number */
  
  for (k=0; k<kmax; k++) {
    x = ix + fx;
    x -= ix; /* careful with compiler optimization */
    if ( ( x != fx ) || ( fabs(ix+fx) >= BIGDBL ) ) break;
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix += 101.0;
    exact = -exact;
  }

  fx = -0.0625;
  exact = 0.980785280403230449126182236134;
  ix = LONG_MIN - 1.0;
  ix -= fabs(fmod(ix,2.0)); /* make sure of even number */
  
  for (k=0; k<kmax; k++) {
    x = ix + fx;
    x -= ix; /* careful with compiler optimization */
    if ( ( x != fx ) || ( fabs(ix+fx) >= BIGDBL ) ) break;
    TEST_SF(s, gsl_sf_cos_pi_e, (ix+fx, &r), exact, TEST_TOL0, GSL_SUCCESS);
    ix -= 101.0;
    exact = -exact;
  }
  
  return s;
}

