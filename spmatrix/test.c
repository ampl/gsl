/* spmatrix/test.c
 * 
 * Copyright (C) 2018 Patrick Alken
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

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spmatrix.h>

int status = 0;

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
#include "test_complex_source.c"
#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
#include "test_complex_source.c"
#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
#include "test_complex_source.c"
#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef  BASE_CHAR

int
main (void)
{
  const size_t num = 3;
  const size_t M[] = { 53, 40, 30 };
  const size_t N[] = { 107, 20, 30 };
  const double density[] = { 0.3, 0.2, 0.5 };
  gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
  size_t i;

  for (i = 0; i < num; ++i)
    {
      test_all (M[i], N[i], density[i], r);
      test_float_all (M[i], N[i], density[i], r);
      test_long_double_all (M[i], N[i], density[i], r);
      test_ulong_all (M[i], N[i], density[i], r);
      test_long_all (M[i], N[i], density[i], r);
      test_uint_all (M[i], N[i], density[i], r);
      test_int_all (M[i], N[i], density[i], r);
      test_ushort_all (M[i], N[i], density[i], r);
      test_short_all (M[i], N[i], density[i], r);
      test_uchar_all (M[i], N[i], density[i], r);
      test_char_all (M[i], N[i], density[i], r);

      test_complex_all (M[i], N[i], density[i], r);
      test_complex_float_all (M[i], N[i], density[i], r);
      test_complex_long_double_all (M[i], N[i], density[i], r);
    }

  gsl_rng_free(r);

  exit (gsl_test_summary ());
}
