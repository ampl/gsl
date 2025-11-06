.. index::
   single: basis splines, B-splines
   single: splines, basis

.. _chap_basis-splines:

*************
Basis Splines
*************

This chapter describes functions for the computation of smoothing
basis splines (B-splines). A smoothing spline differs from an
interpolating spline in that the resulting curve is not required to
pass through each data point.  For information about
interpolating splines, see :ref:`sec_interpolation`.

The header file :file:`gsl_bspline.h` contains the prototypes for the
bspline functions and related declarations.

.. index::
   single: basis splines, overview

Overview
========

Basis splines are commonly used to fit smooth curves to data sets.
They are defined on some interval :math:`[a,b]` which is sub-divided
into smaller intervals by defining a non-decreasing sequence of *knots*.
Between any two adjacent knots, a B-spline of order :math:`k` is a
polynomial of degree :math:`k-1`, and so B-splines defined over the
whole interval :math:`[a,b]` are called *piecewise polynomials*. In
order to generate a smooth spline spanning the full interval, the
different polynomial pieces are often required to match one or more
derivatives at the knots.

B-splines of order :math:`k` are basis functions
for spline functions of the same order, so that all possible spline
functions on some interval :math:`[a,b]` can be expressed as linear combinations
of B-splines,

.. math:: f(x; \mathbf{t}) = \sum_{i=1}^n c_i B_{i,k}(x; \mathbf{t}), \quad a \le x \le b

Here, the :math:`n` coefficients :math:`c_i` are known as
*control points* and are often determined from a least squares fit.
The :math:`B_{i,k}(x; \mathbf{t})` are the basis splines of order :math:`k`,
and depend on a knot vector :math:`\mathbf{t}`, which is a non-decreasing
sequence of :math:`n + k` knots

.. math:: \mathbf{t} = \left\{ t_1, t_2, \dots, t_{n + k} \right\}

On each knot interval :math:`t_i \le x < t_{i+1}`, the B-splines :math:`B_{i,k}(x;\mathbf{t})`
are polynomials of degree :math:`k-1`.
The :math:`n` basis splines of order :math:`k` are defined by

.. only:: not texinfo

   .. math::

      B_{i,1}(x) &=
        \left\{
          \begin{array}{cc}
            1, & t_i \le x < t_{i+1} \\
            0, & \textrm{else}
          \end{array}
        \right. \\
      B_{i,k}(x) &= {(x - t_i) \over (t_{i+k-1} - t_i)} B_{i,k-1}(x) +
                    {(t_{i+k} - x) \over (t_{i+k} - t_{i+1})} B_{i+1,k-1}(x)

.. only:: texinfo

   ::

      B_(i,1)(x) = (1, t_i <= x < t_(i+1)
                   (0, else
      B_(i,k)(x) = [(x - t_i)/(t_(i+k-1) - t_i)] B_(i,k-1)(x)
                    + [(t_(i+k) - x)/(t_(i+k) - t_(i+1))] B_(i+1,k-1)(x)

for :math:`i = 1, \ldots, n`. The common case of cubic B-splines
is given by :math:`k = 4`. The above recurrence relation can be
evaluated in a numerically stable way by the de Boor algorithm.

In what follows, we will drop the :math:`\mathbf{t}` and simply refer
to the spline functions as :math:`f(x)` and :math:`B_{i,k}(x)` for simplicity,
with the understanding that each set of spline basis functions depends on
the chosen knot vector :math:`\mathbf{t}`. We will also occasionally
drop the order :math:`k` and use :math:`B_i(x)` and :math:`B_{i,k}(x)`
interchangeably.

Vector Splines
--------------

The spline above can be generalized to vector control points,

.. math:: \mathbf{f}(x) = \sum_{i=1}^n \mathbf{c}_i B_{i,k}(x)

where the :math:`\mathbf{c}_i` may have some length :math:`m`.

.. index::
   single: basis splines, breakpoints

Knot vector
-----------

In practice, one will often identify locations inside the interval
:math:`[a,b]` where it makes sense to end one polynomial piece and
begin another. These locations are called *interior knots* or *breakpoints*.
It may be desirable to increase the density of knots in regions with rapid
changes, and decrease the density when the underlying data is changing slowly.
For some applications it may be simply best to distribute the knots with
equal spacing on the interval, which are called *uniform knots*.

The full knot vector will look like

.. math:: \mathbf{t} = \{ \underbrace{t_1, \dots, t_{k-1}}_{k-1}, \underbrace{t_k, \dots, t_{n+1},}_{n-k+2 \textrm{ interior knots}} \underbrace{t_{n+2}, \dots, t_{n+k}}_{k-1} \}

The :math:`n-k+2` interior knots are in :math:`[a,b]`. In order to have a full
basis of B-splines, we require an additional :math:`2(k-1)` *exterior knots*,
:math:`\{ t_1,\dots,t_{k-1}\}` and :math:`\{ t_{n+2},\dots,t_{n+k}\}`.
These knots must satisfy

.. math::

   t_1 & \le \cdots \le t_{k-1} = a \\
   b &= t_{n+2} \le \cdots \le t_{n+k}

but can otherwise be arbitrary.

.. index::
   single: basis splines, allocating

Allocation of B-splines
=======================

.. type:: gsl_bspline_workspace

   The computation of B-spline functions requires a preallocated workspace.

.. function:: gsl_bspline_workspace * gsl_bspline_alloc (const size_t k, const size_t nbreak)
              gsl_bspline_workspace * gsl_bspline_alloc_ncontrol (const size_t k, const size_t ncontrol)

   These functions allocate a workspace for computing B-splines of order
   :data:`k`. The number of breakpoints is given by :data:`nbreak`. Alternatively,
   the number of control points :math:`n` may be specified in :data:`ncontrol`.
   The relationship is :math:`\textrm{nbreak} = n - k + 2`. Cubic B-splines
   are specified by :math:`k = 4`. The size of the workspace is
   :math:`O(2k^2 + 5k + \textrm{nbreak})`.

.. function:: void gsl_bspline_free (gsl_bspline_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. index::
   single: knots, basis splines
   single: basis splines, knots

Initializing the knots vector
=============================

.. function:: int gsl_bspline_init_augment (const gsl_vector * tau, gsl_bspline_workspace * w)

   This function initializes the knot vector :math:`\mathbf{t}` using the
   non-decreasing interior knot vector :data:`tau` of length :math:`n - k + 2` by
   choosing the exterior knots to be equal to the first and last interior knot.
   The total knot vector is

   .. math:: \mathbf{t} = \{ \underbrace{\tau_1, \dots, \tau_1}_{k-1 \textrm{ copies}}, \tau_1, \tau_2, \dots, \tau_{n-k+2}, \underbrace{\tau_{n-k+2}, \dots, \tau_{n-k+2}}_{k-1 \textrm{ copies}} \}

   where :math:`\tau_i` is specified by the :math:`i`-th element of :data:`tau`.

   This function does not verify that the :data:`tau` vector is non-decreasing, so
   the user must ensure this.

.. function:: int gsl_bspline_init_uniform (const double a, const double b, gsl_bspline_workspace * w)

   This function constructs a knot vector using uniformly spaced interior knots on :math:`[a,b]` and
   sets the exterior knots to the interval endpoints.
   The total knot vector is

   .. math:: \mathbf{t} = \{ \underbrace{a, \dots, a}_{k-1 \textrm{ copies}}, \tau_1, \tau_2, \dots, \tau_{n-k+2}, \underbrace{b, \dots, b}_{k-1 \textrm{ copies}} \}

   with
   
   .. math:: \tau_i = a + \frac{i-1}{n-k+1} (b-a), \qquad i = 1, \cdots, n-k+2

.. function:: int gsl_bspline_init_periodic (const double a, const double b, gsl_bspline_workspace * w)

   This function constructs a knot vector suitable for defining a periodic B-spline on :math:`[a,b]`
   with a period :math:`P = b - a`. There are several ways to choose a periodic knot vector.
   Here, we adopt the approach of Dierckx, 1982, and require

   .. math::

      \begin{array}{ll}
        t_i = t_{i + n - k + 1} - P, & i = 1, \dots, k-1 \\
        t_{n-k+1+i} = t_i + P, & i = k+1, \dots, 2k-1
      \end{array}

   We further define a uniform spacing between each consecutive knot pair. The knot values are
   given by

   .. math:: t_i = \frac{i-k}{n-k+1} (b-a) + a, \qquad i = 1, \dots, n+k

.. function:: int gsl_bspline_init_interp (const gsl_vector * x, gsl_bspline_workspace * w)

   This function calculates a knot vector suitable for interpolating data at the
   :math:`n` abscissa values provided in :data:`x`, which must be in ascending order.
   The knot vector for interpolation must satisfy the *Schoenberg-Whitney* conditions,

   .. math:: t_i < x_i < t_{i+k}, \quad i = 1, \dots, n

   This function selects the knots according to

   .. math::

      t_i = \left\{
              \begin{array}{ll}
                x_1 & i=1,\dots,k \\
                \frac{1}{k-1} \sum_{j=i-k+1}^{i-1} x_j & i=k+1,\dots,n \\
                x_n & i=n+1,\dots,n+k
              \end{array}
            \right.

   which satisfy the Schoenberg-Whitney conditions, although this may not be
   the optimal knot vector for a given interpolation problem.

.. function:: int gsl_bspline_init_hermite(const size_t nderiv, const gsl_vector * x, gsl_bspline_workspace * w)

   This function calculates a knot vector suitable for performing Hermite
   interpolation data at the :math:`n` abscissa values provided in :data:`x`,
   which must be in ascending order. The parameter :data:`nderiv` specifies
   the maximum derivative order which will be interpolated.

   The knot vector has length :math:`(m+1)(n+2)`, where :math:`m = nderiv`, and is given by,

   .. math:: \mathbf{t} = \{ \underbrace{x_0, \dots, x_0}_{m+1 \textrm{ copies}}, \underbrace{x_1, \dots, x_1}_{m+1}, \dots, \underbrace{x_n, \dots, x_n}_{m+1}, \underbrace{x_{n+1}, \dots, x_{n+1}}_{m+1} \}

   This function uses the convention that :math:`x_0 = x_1` and :math:`x_{n+1} = x_n`. See
   Mummy (1989) for more details.

.. function:: int gsl_bspline_init (const gsl_vector * t, gsl_bspline_workspace * w)

   This function sets the knot vector equal to the user-supplied vector :data:`t` which
   must have length :math:`n+k` and be non-decreasing.

.. index::
   single: basis splines, properties

Properties of B-splines
=======================

.. function:: size_t gsl_bspline_ncontrol (const gsl_bspline_workspace * w)

   This function returns the number of B-spline control points given by
   :math:`n = \textrm{nbreak} + k - 2`.

.. function:: size_t gsl_bspline_order (const gsl_bspline_workspace * w)

   This function returns the B-spline order :math:`k`.

.. function:: size_t gsl_bspline_nbreak (const gsl_bspline_workspace * w)

   This function returns the number of breakpoints (interior knots) defined for the B-spline.

.. index::
   single: basis splines, evaluation

Evaluation of B-splines
=======================

.. function:: int gsl_bspline_calc (const double x, const gsl_vector * c, double * result, gsl_bspline_workspace * w)

   This function evaluates the B-spline

     .. math:: f(x) = \sum_{i=1}^n c_i B_i(x)

   at the point :data:`x` given the :math:`n` control points :data:`c`. The
   result :math:`f(x)` is stored in :data:`result`.
   If :data:`x` lies outside the knot interval, a Taylor series approximation about the closest
   knot endpoint is used to extrapolate the spline.

.. function:: int gsl_bspline_vector_calc (const double x, const gsl_matrix * c, gsl_vector * result, gsl_bspline_workspace * w)

   This function evaluates the vector valued B-spline

     .. math:: \mathbf{f}(x) = \sum_{i=1}^n \mathbf{c}_i B_i(x)

   at the point :data:`x` given the :math:`n` control vectors :math:`\mathbf{c}_i`.
   The vector :math:`\mathbf{c}_i` is stored in the :math:`i`-th column of :data:`c`.
   The result :math:`\mathbf{f}(x)` is stored in :data:`result`. The parameters
   :data:`c` and :data:`result` must have the same number of rows.
   If :data:`x` lies outside the knot interval, a Taylor series approximation about the closest
   knot endpoint is used to extrapolate the spline.

.. function:: int gsl_bspline_basis (const double x, gsl_vector * B, size_t * istart, gsl_bspline_workspace * w)

   This function evaluates all potentially nonzero B-spline basis
   functions at the position :data:`x` and stores them in the vector :data:`B` of length
   :math:`k`, so that the :math:`i`-th element is :math:`B_{(istart+i)}(x)` for
   :math:`i=0,1,\dots,k-1`. The output :data:`istart` is guaranteed to be in
   :math:`[0, n - k]`.
   By returning only the nonzero basis functions, this function
   allows quantities involving linear combinations of the :math:`B_i(x)`
   to be computed without unnecessary terms
   (such linear combinations occur, for example,
   when evaluating an interpolated function).

.. index::
   single: basis splines, derivatives

Evaluation of B-spline derivatives
==================================

.. function:: int gsl_bspline_calc_deriv(const double x, const gsl_vector * c, const size_t nderiv, double * result, gsl_bspline_workspace * w)

   This function evaluates the B-spline derivative

     .. math:: \frac{d^j}{dx^j} f(x) = \sum_{i=1}^n c_i \frac{d^j}{dx^j} B_i(x)

   at the point :data:`x` given the :math:`n` control points :data:`c`. The derivative
   order :math:`j` is specified by :data:`nderiv`. The result
   :math:`d^j f(x) / dx^j` is stored in :data:`result`.
   If :data:`x` lies outside the knot interval, a Taylor series approximation about the closest
   knot endpoint is used to extrapolate the spline.

.. function:: int gsl_bspline_vector_calc_deriv (const double x, const gsl_matrix * c, const size_t nderiv, gsl_vector * result, gsl_bspline_workspace * w)

   This function evaluates the vector valued B-spline derivative

     .. math:: \frac{d^j}{dx^j} \mathbf{f}(x) = \sum_{i=1}^n \mathbf{c}_i \frac{d^j}{dx^j} B_i(x)

   at the point :data:`x` given the :math:`n` control vectors :math:`\mathbf{c}_i`.
   The vector :math:`\mathbf{c}_i` is stored in the :math:`i`-th column of :data:`c`.
   The result :math:`d^j \mathbf{f}(x) / dx^j` is stored in :data:`result`. The parameters
   :data:`c` and :data:`result` must have the same number of rows. The derivative
   order :math:`j` is specified in :data:`nderiv`.
   If :data:`x` lies outside the knot interval, a Taylor series approximation about the closest
   knot endpoint is used to extrapolate the spline.

.. function:: int gsl_bspline_basis_deriv (const double x, const size_t nderiv, gsl_matrix * dB, size_t * istart, gsl_bspline_workspace * w)

   This function evaluates all potentially nonzero B-spline basis function
   derivatives of orders :math:`0` through :data:`nderiv` (inclusive) at
   the position :data:`x` and stores them in the matrix :data:`dB`.  The
   :math:`(i,j)`-th element of :data:`dB` is :math:`d^jB_{(istart+i)}(x)/dx^j`.
   The last row of :data:`dB` contains :math:`d^jB_{istart+k-1}(x)/dx^j`.
   The matrix :data:`dB` must be
   of size :math:`k` by at least :math:`nderiv + 1`.
   The output :data:`istart` is guaranteed to be in :math:`[0, n - k]`.
   Note that function evaluations are included as the zeroth order derivatives in :data:`dB`.
   By returning only the nonzero basis functions, this function allows
   quantities involving linear combinations of the :math:`B_i(x)` and
   their derivatives to be computed without unnecessary terms.

.. index::
   single: basis splines, integration

Evaluation of B-spline integrals
================================

.. function:: int gsl_bspline_calc_integ (const double a, const double b, const gsl_vector * c, double * result, gsl_bspline_workspace * w)

   This function computes the integral

   .. math:: \int_a^b f(x) dx = \sum_{i=1}^n c_i \int_a^b B_i(x) dx

   using the coefficients stored in :data:`c` and the integration limits
   (:data:`a`, :data:`b`). The integral value is stored in :data:`result`
   on output.

.. function:: int gsl_bspline_basis_integ (const double a, const double b, gsl_vector * y, gsl_bspline_workspace * w)

   This function computes the integral from :data:`a` to :data:`b` of each
   basis function :math:`B_i(x)`,

   .. math:: y_i = \int_a^b B_i(x) dx

   storing the results in the vector :data:`y`, of length :code:`ncontrol`.
   The algorithm uses Gauss-Legendre quadrature on each knot interval.

.. index::
   single: basis splines, least squares

Least Squares Fitting with B-splines
====================================

A common application is to fit a B-spline to a set of :math:`m` data
points :math:`(x_i,y_i)` without requiring the spline to pass through
each point, but rather the spline fits the data in a least square sense.
In this case, the number of control points can be much less than :math:`m`
and the result is a smooth spline which minimizes the cost function,

.. math:: \chi^2 = \sum_{i=1}^m w_i \left( y_i - f(x_i) \right)^2 = || y - X c ||_W^2

where :math:`f(x) = \sum_{j=1}^n c_j B_{j,k}(x)` is the B-spline defined above,
:math:`n` is the number of control points for the spline, :math:`W = \textrm{diag}(w)`, and
:math:`X_{ij} = B_{j,k}(x_i)` is the least squares design matrix. The parameters
:math:`w_i` are optional weights which can be assigned to
each data point :math:`(x_i,y_i)`. Because each basis spline :math:`B_{j,k}(x)`
has local support :math:`(B_{j,k}(x) = 0` for :math:`x \notin \left[ t_j,\dots,t_{j+k} \right])`,
the least squares design matrix :math:`X` has only :math:`k`
nonzero entries per row. Therefore, the normal equations matrix :math:`X^T W X` is
:ref:`symmetric and banded <sec_symmetric-banded>`, with lower bandwidth :math:`k - 1`.
The GSL routines below use a banded Cholesky factorization to solve the normal
equations system for the unknown coefficient vector,

.. math:: c = \left( X^T W X \right)^{-1} X^T W y

The normal equations approach is often avoided in least squares problems, since
if :math:`X` is badly conditioned, then :math:`X^T W X` is much worse conditioned,
which can lead to loss of accuracy in the solution vector. However, a theorem
of de Boor [de Boor (2001), XI(4)] states that the condition number of the B-spline basis is independent
of the knot vector :math:`\mathbf{t}` and bounded, depending only on :math:`k`.
Numerical evidence has shown that for reasonably small :math:`k` the normal
equations matrix is well conditioned, and a banded Cholesky factorization can
be used safely.

.. function:: int gsl_bspline_lssolve (const gsl_vector * x, const gsl_vector * y, gsl_vector * c, double * chisq, gsl_bspline_workspace * w)
              int gsl_bspline_wlssolve (const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts, gsl_vector * c, double * chisq, gsl_bspline_workspace * w)

   These functions fit a B-spline to the input data (:data:`x`, :data:`y`) which must be
   the same length :math:`m`. Optional weights may be given in the :data:`wts` vector. The
   banded normal equations matrix is assembled and factored with a Cholesky factorization,
   and then solved, storing the coefficients in the vector :data:`c`, which has length
   :math:`n`. The cost function :math:`\chi^2` is output in :data:`chisq`. If the normal equations matrix
   is singular, the error code :macro:`GSL_EDOM` is returned. This could occur, for example,
   if there is a large gap in the data vectors so that no data points constrain a particular
   :math:`B_{i,k}(x)` basis function.

   On output, :code:`w->XTX` contains the Cholesky factor of the symmetric banded normal equations
   matrix, and this can be passed directly to the functions :func:`gsl_bspline_covariance` and
   :func:`gsl_bspline_rcond`.

.. function:: int gsl_bspline_lsnormal (const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts, gsl_vector * XTy, gsl_matrix * XTX, gsl_bspline_workspace * w)

   This function builds the normal equations matrix :math:`X^T W X` and right hand side vector
   :math:`X^T W y`, which are stored in :data:`XTX` and :data:`XTy` on output. The diagonal
   weight matrix :math:`W` is specified by the input vector :data:`wts`, which may be set
   to :code:`NULL`, in which case :math:`W = I`. The inputs :data:`x`, :data:`y`, and :data:`wts`
   represent the data and weights to be fitted in a least-squares sense, and they all have length
   :math:`m`. While the full normal equations matrix has
   dimension :math:`n`-by-:math:`n`, it has a significant number of zeros, and so the output
   matrix :data:`XTX` is stored in :ref:`banded symmetric format <sec_symmetric-banded>`. Therefore,
   its dimensions are :math:`n`-by-:math:`k`. The vector :data:`XTy` has length :math:`n`.
   The output matrix :data:`XTX` can be passed directly into a banded Cholesky factorization
   routine, such as :func:`gsl_linalg_cholesky_band_decomp`.

   This function is normally called by :func:`gsl_bspline_wlssolve` but is available for
   users with more specialized applications.

.. function:: int gsl_bspline_lsnormalm (const gsl_vector * x, const gsl_matrix * Y, const gsl_vector * wts, gsl_matrix * XTY, gsl_matrix * XTX, gsl_bspline_workspace * w)

   This function builds the normal equations matrix :math:`X^T W X` and right hand side vectors
   :math:`X^T W Y`, which are stored in :data:`XTX` and :data:`XTY` on output. The diagonal
   weight matrix :math:`W` is specified by the input vector :data:`wts`, which may be set
   to :code:`NULL`, in which case :math:`W = I`. The input :data:`Y` is a matrix of the
   right hand side vectors, of size :math:`m`-by-nrhs. The inputs :data:`x` and :data:`wts`
   represent the data and weights to be fitted in a least-squares sense, and they have length
   :math:`m`. While the full normal equations matrix has
   dimension :math:`n`-by-:math:`n`, it has a significant number of zeros, and so the output
   matrix :data:`XTX` is stored in :ref:`banded symmetric format <sec_symmetric-banded>`. Therefore,
   its dimensions are :math:`n`-by-:math:`k`. The matrix :data:`XTY` has size :math:`n`-by-nrhs.
   The output matrix :data:`XTX` can be passed directly into a banded Cholesky factorization
   routine, such as :func:`gsl_linalg_cholesky_band_decomp`.

.. function:: int gsl_bspline_residuals (const gsl_vector * x, const gsl_vector * y, const gsl_vector * c, gsl_vector * r, gsl_bspline_workspace * w)

   This function computes the residual vector :math:`r_i = y_i - f(x_i)` by
   evaluating the spline :math:`f(x)` using the coefficients :data:`c`.

.. function:: int gsl_bspline_covariance (const gsl_matrix * cholesky, gsl_matrix * cov, gsl_bspline_workspace * w)

   This function computes the :math:`n`-by-:math:`n` covariance matrix
   :math:`\left( X^T W X \right)^{-1}` using the banded symmetric Cholesky decomposition
   stored in :data:`cholesky`, which has dimensions :math:`n`-by-:math:`k`.
   The output is stored in :data:`cov`. The function :func:`gsl_linalg_cholesky_band_decomp`
   must be called on the normal equations matrix to compute the Cholesky factor
   prior to calling this function.

.. function:: int gsl_bspline_err (const double x, const size_t nderiv, const gsl_matrix * cov, double * err, gsl_bspline_workspace * w)

   This function computes the standard deviation error of a spline (or its derivative),

   .. math:: f^{(q)}(x) = \sum_{j=1}^n c_j B^{(q)}_j(x)

   at the point :data:`x` using the covariance matrix :data:`cov` calculated previously with
   :func:`gsl_bspline_covariance`. The derivative order is specified by :data:`nderiv`. The
   error is stored in :data:`err` on output, which is defined as,

   .. math:: \delta f = \sqrt{\mathbf{B}^{(q)}(x)^T C \mathbf{B}^{(q)}(x)}

   where :math:`\mathbf{B}^{(q)}(x) = \left( B_1^{(q)}(x), \dots, B_n^{(q)}(x) \right)` is
   the vector of B-spline values evaluated at :math:`x`, and :math:`C` is the covariance matrix.

.. function:: int gsl_bspline_rcond (const gsl_matrix * XTX, double * rcond, gsl_bspline_workspace * w)

   This function estimates the reciprocal condition number (using the 1-norm) of the normal
   equations matrix :math:`X^T W X`, using its Cholesky decomposition, which must be
   computed previously by calling :func:`gsl_linalg_cholesky_band_decomp`.
   The reciprocal condition number estimate, defined as
   :math:`1 / (||X^T W X||_1 \cdot ||(X^T W X)^{-1}||_1)`, is stored in :data:`rcond` on output.

Periodic Spline Fitting
-----------------------

Some applications require fitting a spline to a dataset which has an underlying
periodicity associated with it. In these cases it is desirable to construct a spline
which exhibits the same periodicity. A spline :math:`f(x)` is called periodic on
:math:`[a,b]` if it satisfies the conditions

.. math:: f^{(j)}(a) = f^{(j)}(b), \qquad j = 0, \dots, k - 2

These conditions can be satisfied with the following constraints on the knot
vector and the control points (Dierckx, 1982),

.. math::

   t_i &= t_{i + n - k + 1} - P, \qquad i = 1, \dots, k-1 \\
   t_{n-k+1+i} &= t_i + P, \qquad i = k+1, \dots, 2k-1 \\
   c_i &= c_{n - k + 1 + i}, \qquad i = 1, \dots, k - 1

where :math:`P = b - a` is the period. The constraint on the control
points :math:`c_i` means that there are only :math:`n-k+1` independent
control points, and then the remaining :math:`k-1` are determined. The
independent control points :math:`c^{*} = \{ c_1, \dots, c_{n-k+1} \}`
can be determined by solving the following least squares problem (Dierckx, 1982):

.. math:: \min_{c^{*}} \left|\left| y - X^{*} c^{*} \right|\right|_W^2

where :math:`X^{*}_{ij} = B^{*}_{j,k}(x_i)` and

.. math::

   B^{*}_{j,k}(x) = \left\{
     \begin{array}{ll}
       B_{j,k}(x) + B_{n-k+j+1,k}(x), & j = 1, \dots, k-1 \\
       B_{j,k}(x), & j = k, \dots, n-k+1
     \end{array}
   \right.

.. function:: int gsl_bspline_plssolve (const gsl_vector * x, const gsl_vector * y, gsl_vector * c, double * chisq, gsl_bspline_workspace * w)
              int gsl_bspline_pwlssolve (const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts, gsl_vector * c, double * chisq, gsl_bspline_workspace * w)

   These functions fit a periodic B-spline to the input data (:data:`x`, :data:`y`) which
   must be the same length :math:`m`.  Optional weights may be given in the :data:`wts` vector.
   The input vector :data:`x` must be sorted in ascending order, in order to efficiently
   compute the :math:`QR` decomposition of the least squares system.
   The algorithm computes the :math:`QR` decomposition of the least squares matrix,
   one row at a time using Givens transformations, taking advantage of sparse structure
   due to the local properties of the B-spline basis. On output, the coefficients are stored
   in the vector :data:`c`, which has length :math:`n`. The cost function
   :math:`\chi^2 = || y - X c ||_W^2` is output in :data:`chisq`.

.. function:: int gsl_bspline_plsqr(const gsl_vector * x, const gsl_vector * y, const gsl_vector * wts, gsl_matrix * R, gsl_vector * QTy, double * rnorm, gsl_bspline_workspace * w)

   This function calculates the :math:`R` factor and :math:`Q^T y` vector for a periodic
   spline fitting problem on the input data (:data:`x`, :data:`y`) with optional weights
   :data:`wts`. The input vector :data:`x` must be sorted in ascending order.
   On output, the :math:`(n-k+1)`-by-:math:`(n-k+1)` :math:`R` factor is stored
   in :data:`R`, while the first :math:`n-k+1` elements of the vector :math:`Q^T y`
   are stored in :data:`QTy`. The residual norm :math:`||y - X^{*} c^{*}||_W` is stored in :data:`rnorm`.

.. _sec_bspline-gram:

.. index::
   single: basis splines, Gram matrix

Gram Matrix
===========

Let :math:`\mathbf{B}(x) = (B_1(x), B_2(x), \dots, B_n(x))^T` be the :math:`n`-vector whose
elements are the B-spline basis functions evaluated at :math:`x`. Then
the :math:`n`-by-:math:`n` *outer product matrix* of the B-spline basis at :math:`x`
is defined as

.. math:: A(x) = \mathbf{B}(x) \otimes \mathbf{B}(x) = \mathbf{B}(x) \mathbf{B}(x)^T, \quad A_{ij}(x) = B_i(x) B_j(x)

This can be generalized to higher order derivatives,

.. math:: A_{ij}^{(q)}(x) = \left( \frac{d^q}{dx^q} B_i(x) \right) \left( \frac{d^q}{dx^q} B_j(x) \right)

Then, the :math:`n`-by-:math:`n` Gram matrix of the B-spline basis is defined as

.. math:: G_{ij}^{(q)} = \int_a^b A_{ij}^{(q)}(x) dx = \int_a^b \left( \frac{d^q}{dx^q} B_i(x) \right) \left( \frac{d^q}{dx^q} B_j(x) \right) dx

where the integral is usually taken over the full support of the B-splines, i.e. from
:math:`\min(\mathbf{t})` to :math:`\max(\mathbf{t})`. However it can be computed
over any interval :math:`[a,b]`.

The B-spline outer product and Gram matrices are :ref:`symmetric and banded <sec_symmetric-banded>` with lower bandwidth
:math:`k - 1`. These matrices are positive semi-definite, and in the case of :math:`q = 0`,
:math:`A^{(0)}` and :math:`G^{(0)}` are positive
definite. These matrices have applications in least-squares fitting.

.. function:: int gsl_bspline_oprod (const size_t nderiv, const double x, gsl_matrix * A, gsl_bspline_workspace * w)

   This function calculates the outer product matrix :math:`A^{(q)}(x)` at the location :data:`x`, where the
   derivative order :math:`q` is specified by :data:`nderiv`. The matrix is stored in :data:`A` in
   :ref:`symmetric banded <sec_symmetric-banded>` format, and so :data:`A` has dimensions
   :math:`n`-by-:math:`k`. Since the B-spline functions have local support, many elements of the
   outer product matrix could be zero.

.. function:: int gsl_bspline_gram (const size_t nderiv, gsl_matrix * G, gsl_bspline_workspace * w)

   This function calculates the positive semi-definite Gram matrix :math:`G_{ij}^{(q)}`, where the
   derivative order :math:`q` is specified by :data:`nderiv`. The matrix is stored in :data:`G` in
   :ref:`symmetric banded <sec_symmetric-banded>` format, and so :data:`G` has dimensions
   :math:`n`-by-:math:`k`. The integrals are computed over the full support of the B-spline basis.

   The algorithm uses Gauss-Legendre quadrature on each knot interval. See algorithm 5.22 of
   L. Schumaker, "Spline Functions: Basic Theory".

.. function:: int gsl_bspline_gram_interval (const double a, const double b, const size_t nderiv, gsl_matrix * G, gsl_bspline_workspace * w)

   This function calculates the positive semi-definite Gram matrix :math:`G_{ij}^{(q)}`, where the
   derivative order :math:`q` is specified by :data:`nderiv`. The matrix is stored in :data:`G` in
   :ref:`symmetric banded <sec_symmetric-banded>` format, and so :data:`G` has dimensions
   :math:`n`-by-:math:`k`. The integrals are computed over the interval :math:`[a,b]`.

   The algorithm uses Gauss-Legendre quadrature on each knot interval. See algorithm 5.22 of
   L. Schumaker, "Spline Functions: Basic Theory".

.. index::
   single: basis splines, interpolation

Interpolation with B-splines
============================

Given a set of :math:`n` data points, :math:`(x_i,y_i)`, we can
define an interpolating spline with :math:`n` control points
which passes through each data point,

.. math:: f(x_i) = \sum_{j=1}^n c_j B_j(x_i) = y_i, \quad i = 1,\dots,n

This is a square system of equations which can be solved by inverting
the :math:`n`-by-:math:`n` *collocation matrix*, whose elements are
:math:`B_j(x_i)`. In order for the collocation matrix to be invertible,
the :math:`x_i` and spline knots must satisfy the Schoenberg-Whitney conditions,
which are

.. math:: t_i < x_i < t_{i+k}, \quad i = 1,\dots,n

If these conditions are met, then the collocation matrix is banded,
with both the upper and lower bandwidths equal to :math:`k-1`. This
property leads to savings in storage and computation when solving
interpolation problems. To construct a knot vector which satisfies
these conditions, see the function :func:`gsl_bspline_init_interp`.

.. function:: int gsl_bspline_col_interp(const gsl_vector * x, gsl_matrix * XB, gsl_bspline_workspace * w)

   This function constructs the banded collocation matrix to interpolate
   the :math:`n` abscissa values provided in :data:`x`. The values in
   :data:`x` must be sorted in ascending order, although this is not
   checked explicitly by the function. The collocation
   matrix is stored in :data:`XB` on output in banded storage format suitable
   for use with the banded LU routines. Therefore, :data:`XB` must have dimension
   :math:`n`-by-:math:`3(k-1) + 1`.

.. index::
   single: basis splines, Hermite interpolation

Hermite Interpolation with B-splines
====================================

Given a set of :math:`n` data interpolation sites,
:math:`x_1 < x_2 < \dots < x_n`, and a set of interpolation
values,

.. math:: f_i^{(j)}, \qquad i=1,\dots,n, j=0,1,\dots,m

then there is a unique interpolating spline of order
:math:`k = 2m+2` which will match the given function values
(and their derivatives). The spline has :math:`(m+1)n` control
points and is given by,

.. math:: f(x) = \sum_{j=1}^{(m+1)n} c_j B_j(x)

This spline will satisfy,

.. math:: f^{(j)}(x_i) = f_i^{(j)}, \qquad i=1,\dots,n, j=0,1,\dots,m

Mummy (1989) provides an efficient algorithm to determine the spline
coefficients :math:`c_j`.

.. function:: int gsl_bspline_interp_chermite(const gsl_vector * x, const gsl_vector * y, const gsl_vector * dy, gsl_vector * c, const gsl_bspline_workspace * w)

   This function will calculate the coefficients of an interpolating
   cubic Hermite spline which match the provided function values :math:`y_i`
   and first derivatives :math:`y'_i` at the interpolation sites
   :math:`x_i`. The inputs :data:`x`, :data:`y`, :data:`dy` are
   all of the same length :math:`n` and contain the :math:`x_i`,
   :math:`y_i`, and :math:`y'_i` respectively. On output, the spline coefficients
   are stored in :data:`c`, which must have length :math:`2n`.

   This function is specifically designed for cubic Hermite splines, and
   so the spline order must be :math:`k = 4`. The function
   :func:`gsl_bspline_init_hermite` must be called to initialize the knot vector.
   The algorithm used is provided by Mummy (1989).

Projection onto the B-spline Basis
==================================

A function :math:`f(x)` can be projected onto the subspace :math:`V = \textrm{span}(B_i)` spanned
by the B-splines as follows,

.. only:: not texinfo

    .. math::

       f(x) \approx \textrm{proj}_V(f) = \sum_j c_j B_j(x)

.. only:: texinfo

   ::
   
      f(x) =~ proj_V(f) = \sum_j c_j B_j(x)

If :math:`f(x)` is a piecewise polynomial of degree less than :math:`k`, then :math:`f(x)` will
lie in :math:`V` and :math:`f(x) = \textrm{proj}_V(f)`. The expansion coefficients :math:`c_j` can be
calculated as follows,

.. only:: not texinfo

   .. math::

      f(x) B_i(x) &= \sum_j c_j B_j(x) B_i(x) \\
      \underbrace{\int f(x) B_i(x) dx}_{y_i} &= \sum_j c_j \underbrace{\int B_i(x) B_j(x) dx}_{G_{ij}} \\
      G c &= y

.. only:: texinfo

   ::

      f(x) B_i(x) = \sum_j c_j B_j(x) B_i(x)
      \int f(x) B_i(x) dx = \sum_j c_j \int B_i(x) B_j(x) dx
      G c = y

where :math:`G` is the :ref:`B-spline Gram Matrix <sec_bspline-gram>`. The above equation can
be solved using a Cholesky factorization of :math:`G`.

.. function:: int gsl_bspline_proj_rhs (const gsl_function * F, gsl_vector * y, gsl_bspline_workspace * w)

   This function computes the right hand side vector :math:`y_i = \int f(x) B_i(x) dx` for projecting
   a function :data:`F` onto the B-spline basis. The output is stored in :data:`y` which must have
   size :code:`ncontrol`. The integrals are computed using Gauss-Legendre quadrature, which will
   produce exact results if :math:`f(x)` is a piecewise polynomial of degree less than :math:`k`.

.. index::
   single: basis splines, Greville abscissae
   single: basis splines, Marsden-Schoenberg points

Greville abscissae
==================

The Greville abscissae are defined to be the mean location of :math:`k-1`
consecutive knots in the knot vector for each basis spline function of order
:math:`k`.  With the first and last knots in the knot vector excluded,
there are :math:`n` Greville abscissae
for any given B-spline basis.  These values are often used in B-spline
collocation applications and may also be called Marsden-Schoenberg points.

.. function:: double gsl_bspline_greville_abscissa (const size_t i, gsl_bspline_workspace * w)

   Returns the location of the :math:`i`-th Greville abscissa for the given
   B-spline basis.  For the ill-defined case when :math:`k = 1`, the implementation
   chooses to return breakpoint interval midpoints.

.. See https://savannah.gnu.org/bugs/index.php?34361
.. @deftypefun int gsl_bspline_knots_greville (const gsl_vector * abscissae, gsl_bspline_workspace * w, double * abserr);
.. Given target Greville abscissae values in :data:`abscissae` and a workspace
.. :data:`w` where @code{abscissae->size == gsl_bspline_ncontrol(w)}, this functions
.. computes and stores the knots required for the workspace to best approximate
.. the target abscissae.  The approximation is optimal in that the first and last
.. values in :data:`abscissae` are preserved exactly while the 2-norm of the error
.. in any other abscissae is minimized.  If not-@code{NULL}, the sum of the
.. absolute approximation errors over each abscissa is returned in :data:`abserr`.
..
.. The workspace order must satisfy :math:`k > 1` and :data:`abscissae` should be
.. monotonically increasing.  Beware that when @code{w->nbreak} is small relative
.. to @code{w->k} the best approximation may still be of poor quality for
.. non-uniformly spaced :data:`abscissae`.  This function has memory and runtime
.. overhead that scales like a QR-based linear least squares solution on a
.. @code{(abscissae->size - 2)} by @code{(w->nbreak - 2)} problem.
.. @end deftypefun

.. index::
   single: basis splines, examples

Examples
========

Example 1: Basis Splines and Uniform Knots
------------------------------------------

The following program computes and outputs linear, quadratic, cubic, and
quartic basis splines using uniform knots on the interval :math:`[0,1]`.
The knot vector is

.. math:: \mathbf{t} = \left( \underbrace{0, \dots, 0}_{k}, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, \underbrace{1, \dots, 1}_{k} \right)

where the first and last knot are repeated :math:`k` times, where :math:`k` is the spline
order.  The output is shown in :numref:`fig_bspline_knots`.

.. _fig_bspline_knots:

.. figure:: /images/bspline_knots.png
   :scale: 60%

   Linear, quadratic, cubic, and quartic basis splines defined on the interval
   :math:`[0,1]` with uniform knots. Black circles denote the location of
   the knots.

This figure makes it clear that each B-spline has local support and spans
:math:`k + 1` knots. The B-splines near the interval endpoints appear to span
less knots, but note that the endpoint knots have multiplicities of :math:`k`.

The program is given below.

.. include:: examples/bspline_knots.c
   :code:

Example 2: Derivatives of Basis Splines
---------------------------------------

The following program computes and outputs cubic B-splines and their
derivatives using :math:`6` breakpoints and uniform knots on the interval
:math:`[0,1]`. All derivatives up to order :math:`3` are computed.
The output is shown in :numref:`fig_bspline2`.

.. _fig_bspline2:

.. figure:: /images/bspline_deriv.png
   :scale: 60%

   Cubic B-splines and derivatives up to order 3.  Black circles
   denote the location of the knots.

The program is given below.

.. include:: examples/bspline_deriv.c
   :code:

Example 3: Least squares and breakpoints
----------------------------------------

The following program computes a linear least squares fit to data using
cubic B-splines with uniform breakpoints. The data is generated from
the curve :math:`y(x) = \cos{x} \exp{\left( -x/10 \right)}` with Gaussian noise added.
Splines are fitted with :math:`10` and :math:`40` breakpoints for comparison.

The program output is shown below::

  $ ./a.out > bspline_lsbreak.txt
  40 breakpoints: chisq/dof = 1.008999e+00
  10 breakpoints: chisq/dof = 1.014761e+00

The data and fitted models are shown in :numref:`fig_bspline_lsbreak`.

.. _fig_bspline_lsbreak:

.. figure:: /images/bspline_lsbreak.png
   :scale: 60%

   Cubic B-spline least squares fit

The program is given below.

.. include:: examples/bspline_lsbreak.c
   :code:

Example 4: Least squares and spline order
-----------------------------------------

The following program computes least squares fits to the same
dataset in the previous example, using :math:`10` uniform
breakpoints and computing splines of order :math:`1,2,3,4,5`.

The data and fitted models are shown in :numref:`fig_bspline_lsorder`.

.. _fig_bspline_lsorder:

.. figure:: /images/bspline_lsorder.png
   :scale: 60%

   B-spline least squares fit with varying spline orders

The program is given below.

.. include:: examples/bspline_lsorder.c
   :code:

Example 5: Least squares, regularization and the Gram matrix
------------------------------------------------------------

The following example program fits a cubic B-spline to data generated
from the Gaussian

.. math:: g(x) = e^{-x^2}

on the interval :math:`[-1.5,1.5]` with noise added. Two data gaps are
deliberately introduced on :math:`[-1.1,-0.7]` and :math:`[0.1,0.55]`,
which makes the problem ill-posed. :numref:`fig_bspline_gram1` shows
that the ordinary least squares solution exhibits large errors in
the data gap regions due to the ill-conditioned least squares matrix.

One way to correct the situation is to minimize the data misfit but
also minimize the curvature of the spline averaged over the full interval.
We can define the average curvature using the second derivative of the
spline as follows:

.. math::
   
   \textrm{curvature} &\approx \int \left| \frac{d^2}{dx^2} f(x) \right|^2 dx \\
                      &= \int \left| \frac{d^2}{dx^2} \sum_{i=1}^n c_i B_i(x) \right|^2 dx \\
                      &= \int \left| \sum_{i=1}^n c_i B''_i(x) \right|^2 dx \\
                      &= \int \left| c^T B''(x) \right|^2 dx \\
                      &= \int c^T B''(x) B''^T(x) c dx \\
                      &= c^T \left( \int B''(x) B''^T(x) dx \right) c \\
                      &= c^T G^{(2)} c

where :math:`G^{(2)}` is the :ref:`Gram matrix <sec_bspline-gram>` of second derivatives of the B-splines.
Therefore, our regularized least squares problem is

.. math:: \min_{c} || y - X c ||_W^2 + \lambda^2 c^T G^{(2)} c

where the regularization parameter :math:`\lambda^2` represents a tradeoff between
minimizing the data misfit and minimizing the spline curvature.
The solution of this least-squares problem is

.. math:: c = \left( X^T W X + \lambda^2 G^{(2)} \right)^{-1} X^T W y

The example program below solves this problem without regularization :math:`(\lambda^2 = 0)`
and with regularization :math:`(\lambda^2 = 0.1)` for comparison. We also plot
the :math:`1 \sigma` confidence intervals for both solutions using the diagonal elements
of the covariance matrix :math:`\textrm{Cov} = \left( X^T W X + \lambda^2 G^{(2)} \right)^{-1}`.

:numref:`fig_bspline_gram1` shows the result of regularizing the spline over the full
interval, using the command

::

  ./a.out > data.txt

.. _fig_bspline_gram1:

.. figure:: /images/bspline_gram1.png
   :scale: 60%

   B-spline least squares fit with exact (green), unregularized (red with :math:`1 \sigma` confidence intervals),
   and fully regularized (blue with :math:`1 \sigma` confidence intervals) solutions.

:numref:`fig_bspline_gram2` shows the result of regularizing the spline over the smaller
interval :math:`[0,0.6]`, correcting only the second gap, using the command

::

  ./a.out 0 0.6 > data.txt

.. _fig_bspline_gram2:

.. figure:: /images/bspline_gram2.png
   :scale: 60%

   B-spline least squares fit with exact (green), unregularized (red with :math:`1 \sigma` confidence intervals),
   and partially regularized (blue with :math:`1 \sigma` confidence intervals) solutions.

The program is given below.

.. include:: examples/bspline_gram.c
   :code:

Example 6: Least squares and endpoint regularization
----------------------------------------------------

Data fitting with splines can often produce undesired effects at the spline
endpoints, where they are less constrained by data. As an example, we fit
order :math:`k = 10` splines to the Runge function,

.. math:: g(x) = \frac{1}{1 + 25 x^2}

with 20 uniform breakpoints on the interval :math:`[a,b] = [-1,1]`.
:numref:`fig_bspline_lsend` shows that the fitted spline (red) exhibits
a sharp upward feature at the left endpoint. This is due to the lack
of data constraining this part of the spline, and would likely result
in errors if the user wanted to extrapolate the spline outside the fit
interval. One method to correct this is to apply regularization to
the spline at the endpoints by minimizing one or more of its derivatives.
In this example, we minimize the first derivative of the spline at both
endpoints, :math:`|f'(a)|^2` and :math:`|f'(b)|^2`. Note that

.. math::

   \left| \frac{d}{dx} f(x) \right|^2 &= \left| \frac{d}{dx} \sum_{i=1}^n c_i B_i(x) \right|^2 \\
   &= \left| c^T B'(x) \right|^2 \\
   &= c^T B'(x) B'^T(x) c \\
   &= c^T A^{(1)}(x) c

where :math:`A^{(1)}` is the :ref:`outer product matrix <sec_bspline-gram>` of first derivatives.
Therefore, our regularized least squares problem is

.. math:: \min_{c} || y - X c ||_W^2 + \lambda^2 c^T A^{(1)}(a) c + \lambda^2 c^T A^{(1)}(b) c

where we have chosen to use the same regularization parameter :math:`\lambda^2` to regularize both
ends of the spline. In the example program below, we compute the unregularized solution (:math:`\lambda^2 = 0`)
and regularized solution (:math:`\lambda^2 = 10`). The program outputs the first derivative of
the spline at both endpoints::

  unregularized endpoint deriv 1: [ -1.081170e+01,  -2.963725e+00]
    regularized endpoint deriv 1: [ -2.735857e-02,  -5.371758e-03]

showing that the regularization has forced the first derivative smaller at both endpoints. The
results are shown in :numref:`fig_bspline_lsend`.

.. _fig_bspline_lsend:

.. figure:: /images/bspline_lsend.png
   :scale: 60%

   Runge function B-spline least squares fit with exact (green), unregularized (red) and regularized (blue)
   solutions. The zoomed inset views show the spline endpoints with :math:`1 \sigma` confidence intervals.

The program is given below.

.. include:: examples/bspline_lsend.c
   :code:

Example 7: Least squares and periodic splines
---------------------------------------------

The following example program fits B-spline curves to a dataset which is periodic.
Both a non-periodic and periodic spline are fitted to the data for comparison. The
data is derived from the function on :math:`[0,2\pi]`,

.. math:: f(x) = \sin{x} - \cos{2x}

with Gaussian noise added. The program fits order :math:`k=6` splines to the
dataset.

.. _fig_bspline_per:

.. figure:: /images/bspline_per.png
   :scale: 60%

   Data (black) with non-periodic spline fit (red) and periodic spline fit (green).

The program additionally computes the derivatives up to order :math:`5` at
the endpoints :math:`x=0` and :math:`x=2\pi`::

  === Non-periodic spline endpoint derivatives ===
  deriv 0: [ -8.697939e-01,  -8.328553e-01]
  deriv 1: [ -2.423132e+00,   4.277439e+00]
  deriv 2: [  3.904362e+01,   3.710592e+01]
  deriv 3: [ -2.096142e+02,   2.036125e+02]
  deriv 4: [  7.240121e+02,   7.384888e+02]
  deriv 5: [ -1.230036e+03,   1.294264e+03]
  === Periodic spline endpoint derivatives ===
  deriv 0: [ -1.022798e+00,  -1.022798e+00]
  deriv 1: [  1.041910e+00,   1.041910e+00]
  deriv 2: [  4.466384e+00,   4.466384e+00]
  deriv 3: [ -1.356654e+00,  -1.356654e+00]
  deriv 4: [ -2.751544e+01,  -2.751544e+01]
  deriv 5: [  5.055419e+01,   2.674261e+01]

We can see that the periodic fit matches derivative values up
to order :math:`4`, with the :math:`5`-th derivative being non-continuous.

The program is given below.

.. include:: examples/bspline_per.c
   :code:

Example 8: Projecting onto the B-spline basis
---------------------------------------------

The following example program projects the function

.. math:: f(x) = 3 x^3 - 2 x^2 - 7 x

onto a cubic B-spline basis with uniform knots on the interval :math:`[-2,2]`.
Because the function is cubic, it can be exactly decomposed in the B-spline basis.

.. _fig_bspline_proj:

.. figure:: /images/bspline_proj.png
   :scale: 60%

   Exact function (red) and spline projection (green). Small offset deliberately added to
   green curve.

The program is given below.

.. include:: examples/bspline_proj.c
   :code:

Example 9: Interpolation
------------------------

The following example program computes an interpolating
cubic B-spline for a set of data. It constructs the collocation
matrix for the dataset and factors it with a banded LU decomposition.
The results are shown in :numref:`fig_bspline_interp`.

.. _fig_bspline_interp:

.. figure:: /images/bspline_interp.png
   :scale: 60%

   Data points are shown in black and interpolated spline in light blue.

The program is given below.

.. include:: examples/bspline_interp.c
   :code:

References and Further Reading
==============================

Further information on the algorithms described in this section can be
found in the following references,

* C. de Boor, *A Practical Guide to Splines* (2001), Springer-Verlag,
  ISBN 0-387-95366-3.

* L. L. Schumaker, *Spline Functions: Basic Theory*, 3rd edition,
  Cambridge University Press, 2007.

* P. Dierckx, *Algorithms for smoothing data with periodic and parametric splines*,
  Computer Graphics and Image Processing, 20, 1982.

* P. Dierckx, *Curve and surface fitting with splines*, Oxford University Press, 1995.

* M. S. Mummy, *Hermite interpolation with B-splines*,
  Computer Aided Geometric Design, 6, 177--179, 1989.

Further information of Greville abscissae and B-spline collocation
can be found in the following paper,

* Richard W. Johnson, Higher order B-spline collocation at the Greville
  abscissae.  *Applied Numerical Mathematics*. vol.: 52, 2005, 63--75.

A large collection of B-spline routines is available in the
PPPACK library available at http://www.netlib.org/pppack,
which is also part of SLATEC.
