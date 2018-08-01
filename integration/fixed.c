/* integration/fixed.c
 * 
 * Copyright (C) 2017 Patrick Alken
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

/* the code in this module performs fixed-point quadrature calculations for
 * integrands and is based on IQPACK */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

static int fixed_compute(const double a, const double b, const double alpha, const double beta,
                         gsl_integration_fixed_workspace * w);
static int imtqlx ( const int n, double d[], double e[], double z[] );

gsl_integration_fixed_workspace *
gsl_integration_fixed_alloc(const gsl_integration_fixed_type * type, const size_t n,
                            const double a, const double b, const double alpha, const double beta)
{
  int status;
  gsl_integration_fixed_workspace *w;

  /* check inputs */
  if (n < 1)
    {
      GSL_ERROR_VAL ("workspace size n must be at least 1", GSL_EDOM, 0);
    }

  w = calloc(1, sizeof(gsl_integration_fixed_workspace));
  if (w == NULL)
    {
      GSL_ERROR_VAL ("unable to allocate workspace", GSL_ENOMEM, 0);
    }

  w->weights = malloc(n * sizeof(double));
  if (w->weights == NULL)
    {
      gsl_integration_fixed_free(w);
      GSL_ERROR_VAL ("unable to allocate weights", GSL_ENOMEM, 0);
    }

  w->x = malloc(n * sizeof(double));
  if (w->x == NULL)
    {
      gsl_integration_fixed_free(w);
      GSL_ERROR_VAL ("unable to allocate x", GSL_ENOMEM, 0);
    }

  w->diag = malloc(n * sizeof(double));
  if (w->diag == NULL)
    {
      gsl_integration_fixed_free(w);
      GSL_ERROR_VAL ("unable to allocate diag", GSL_ENOMEM, 0);
    }

  w->subdiag = malloc(n * sizeof(double));
  if (w->subdiag == NULL)
    {
      gsl_integration_fixed_free(w);
      GSL_ERROR_VAL ("unable to allocate subdiag", GSL_ENOMEM, 0);
    }

  w->n = n;
  w->type = type;

  /* compute quadrature weights and nodes */
  status = fixed_compute(a, b, alpha, beta, w);
  if (status)
    {
      gsl_integration_fixed_free(w);
      GSL_ERROR_VAL ("error in integration parameters", GSL_EDOM, 0);
    }

  return w;
}

void
gsl_integration_fixed_free(gsl_integration_fixed_workspace * w)
{
  if (w->weights)
    free(w->weights);

  if (w->x)
    free(w->x);

  if (w->diag)
    free(w->diag);

  if (w->subdiag)
    free(w->subdiag);

  free(w);
}

size_t
gsl_integration_fixed_n(const gsl_integration_fixed_workspace * w)
{
  return w->n;
}

double *
gsl_integration_fixed_nodes(const gsl_integration_fixed_workspace * w)
{
  return w->x;
}

double *
gsl_integration_fixed_weights(const gsl_integration_fixed_workspace * w)
{
  return w->weights;
}

int
gsl_integration_fixed(const gsl_function * func, double * result,
                      const gsl_integration_fixed_workspace * w)
{
  const size_t n = w->n;
  size_t i;
  double sum = 0.0;

  for (i = 0; i < n; ++i)
    {
      double fi = GSL_FN_EVAL(func, w->x[i]);
      sum += w->weights[i] * fi;
    }

  *result = sum;

  return GSL_SUCCESS;
}

/*
fixed_compute()
  Compute quadrature weights and nodes
*/

static int
fixed_compute(const double a, const double b, const double alpha, const double beta,
              gsl_integration_fixed_workspace * w)
{
  int s;
  const size_t n = w->n;
  gsl_integration_fixed_params params;
  size_t i;

  params.a = a;
  params.b = b;
  params.alpha = alpha;
  params.beta = beta;

  /* check input parameters */
  s = (w->type->check)(n, &params);
  if (s)
    return s;

  /* initialize Jacobi matrix */
  s = (w->type->init)(n, w->diag, w->subdiag, &params);
  if (s)
    return s;

  if (params.zemu <= 0.0)
    {
      GSL_ERROR("zeroth moment must be positive", GSL_EINVAL);
    }

  for ( i = 0; i < n; i++ )
    {
      w->x[i] = w->diag[i];
    }

  w->weights[0] = sqrt (params.zemu);

  for ( i = 1; i < n; i++ )
    {
      w->weights[i] = 0.0;
    }

  /* diagonalize the Jacobi matrix */
  s = imtqlx (n, w->x, w->subdiag, w->weights);
  if (s)
    return s;

  for (i = 0; i < n; i++)
    {
      w->weights[i] = w->weights[i] * w->weights[i];
    }

  /*
   * The current weights and nodes are valid for a = 0, b = 1.
   * Now scale them for arbitrary a,b
   */
  {
    double p = pow ( params.slp, params.al + params.be + 1.0 );
    size_t k;

    for ( k = 0; k < n; k++ )
      {
        w->x[k] = params.shft + params.slp * w->x[k];
        w->weights[k] = w->weights[k] * p;
      }
  }

  return GSL_SUCCESS;
}

/******************************************************************************/
/*
  Purpose:

    IMTQLX diagonalizes a symmetric tridiagonal matrix.

  Discussion:

    This routine is a slightly modified version of the EISPACK routine to 
    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 

    The authors thank the authors of EISPACK for permission to use this
    routine. 

    It has been modified to produce the product Q' * Z, where Z is an input 
    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
    The changes consist (essentially) of applying the orthogonal transformations
    directly to Z as they are generated.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

    Roger Martin, James Wilkinson,
    The Implicit QL Algorithm,
    Numerische Mathematik,
    Volume 12, Number 5, December 1968, pages 377-383.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D(N), the diagonal entries of the matrix.
    On output, the information in D has been overwritten.

    Input/output, double E(N), the subdiagonal entries of the 
    matrix, in entries E(1) through E(N-1).  On output, the information in
    E has been overwritten.

    Input/output, double Z(N).  On input, a vector.  On output,
    the value of Q' * Z, where Q is the matrix that diagonalizes the
    input symmetric tridiagonal matrix.
*/

static int
imtqlx ( const int n, double d[], double e[], double z[] )
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double r;
  double s;

  if ( n == 1 )
  {
    return GSL_SUCCESS;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( fabs ( e[m-1] ) <= GSL_DBL_EPSILON * ( fabs ( d[m-1] ) + fabs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        return GSL_EMAXITER;
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + fabs ( r ) * GSL_SIGN ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( fabs ( g ) <= fabs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
/*
  Sorting.
*/
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }

  return GSL_SUCCESS;
}
