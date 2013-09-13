#include <config.h>
#include <math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>

static int compute_rank(gsl_vector *v);

#define BASE_LONG_DOUBLE
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
#include "covariance_source.c"
#include "templates_off.h"
#undef  BASE_CHAR

/*
compute_rank()
  Compute rank of a sorted vector

Inputs: v - sorted data vector on input; rank vector on output

Notes: ranks are always computed in double precision
*/

static int
compute_rank(gsl_vector *v)
{
  const size_t n = v->size;
  size_t i = 0;

  while (i < n - 1)
    {
      double vi = gsl_vector_get(v, i);

      if (vi == gsl_vector_get(v, i + 1))
        {
          size_t j = i + 2;
          size_t k;
          double rank = 0.0;

          /* we have detected a tie, find number of equal elements */
          while (j < n && vi == gsl_vector_get(v, j))
            ++j;

          /* compute rank */
          for (k = i; k < j; ++k)
            rank += k + 1.0;

          /* divide by number of ties */
          rank /= (double) (j - i);

          for (k = i; k < j; ++k)
            gsl_vector_set(v, k, rank);

          i = j;
        }
      else
        {
          /* no tie - set rank to natural ordered position */
          gsl_vector_set(v, i, i + 1.0);
          ++i;
        }
    }

  if (i == n - 1)
    gsl_vector_set(v, n - 1, (double) n);

  return GSL_SUCCESS;
} /* compute_rank() */
