.. index::
   single: minimization, multidimensional

*****************************
Multidimensional Minimization
*****************************

This chapter describes routines for finding minima of arbitrary
multidimensional functions.  The library provides low level components
for a variety of iterative minimizers and convergence tests.  These can
be combined by the user to achieve the desired solution, while providing
full access to the intermediate steps of the algorithms.  Each class of
methods uses the same framework, so that you can switch between
minimizers at runtime without needing to recompile your program.  Each
instance of a minimizer keeps track of its own state, allowing the
minimizers to be used in multi-threaded programs. The minimization
algorithms can be used to maximize a function by inverting its sign.

The header file :file:`gsl_multimin.h` contains prototypes for the
minimization functions and related declarations.  

Overview
========

The problem of multidimensional minimization requires finding a point
:math:`x` such that the scalar function,

.. math:: f(x_1, \dots, x_n)

takes a value which is lower than at any neighboring point. For smooth
functions the gradient :math:`g = \nabla f` vanishes at the minimum. In
general there are no bracketing methods available for the
minimization of :math:`n`-dimensional functions.  The algorithms
proceed from an initial guess using a search algorithm which attempts
to move in a downhill direction. 

Algorithms making use of the gradient of the function perform a
one-dimensional line minimisation along this direction until the lowest
point is found to a suitable tolerance.  The search direction is then
updated with local information from the function and its derivatives,
and the whole process repeated until the true :math:`n`-dimensional
minimum is found.

Algorithms which do not require the gradient of the function use
different strategies.  For example, the Nelder-Mead Simplex algorithm
maintains :math:`n+1` trial parameter vectors as the vertices of a
:math:`n`-dimensional simplex.  On each iteration it tries to improve
the worst vertex of the simplex by geometrical transformations.  The
iterations are continued until the overall size of the simplex has
decreased sufficiently.

Both types of algorithms use a standard framework. The user provides a
high-level driver for the algorithms, and the library provides the
individual functions necessary for each of the steps.  There are three
main phases of the iteration.  The steps are,

* initialize minimizer state, :data:`s`, for algorithm :data:`T`
* update :data:`s` using the iteration :data:`T`
* test :data:`s` for convergence, and repeat iteration if necessary

Each iteration step consists either of an improvement to the
line-minimisation in the current direction or an update to the search
direction itself.  The state for the minimizers is held in a
:type:`gsl_multimin_fdfminimizer` struct or a
:type:`gsl_multimin_fminimizer` struct.

.. index::
   single: Multimin, caveats

Caveats
=======

Note that the minimization algorithms can only search for one local
minimum at a time.  When there are several local minima in the search
area, the first minimum to be found will be returned; however it is
difficult to predict which of the minima this will be.  In most cases,
no error will be reported if you try to find a local minimum in an area
where there is more than one.

It is also important to note that the minimization algorithms find local
minima; there is no way to determine whether a minimum is a global
minimum of the function in question.

Initializing the Multidimensional Minimizer
===========================================

The following function initializes a multidimensional minimizer.  The
minimizer itself depends only on the dimension of the problem and the
algorithm and can be reused for different problems.

.. type:: gsl_multimin_fdfminimizer

   This is a workspace for minimizing functions using derivatives.

.. type:: gsl_multimin_fminimizer

   This is a workspace for minimizing functions without derivatives.

.. function:: gsl_multimin_fdfminimizer * gsl_multimin_fdfminimizer_alloc (const gsl_multimin_fdfminimizer_type * T, size_t n)
              gsl_multimin_fminimizer * gsl_multimin_fminimizer_alloc (const gsl_multimin_fminimizer_type * T, size_t n)

   This function returns a pointer to a newly allocated instance of a
   minimizer of type :data:`T` for an :data:`n`-dimension function.  If there
   is insufficient memory to create the minimizer then the function returns
   a null pointer and the error handler is invoked with an error code of
   :macro:`GSL_ENOMEM`.

.. function:: int gsl_multimin_fdfminimizer_set (gsl_multimin_fdfminimizer * s, gsl_multimin_function_fdf * fdf, const gsl_vector * x, double step_size, double tol)
              int gsl_multimin_fminimizer_set (gsl_multimin_fminimizer * s, gsl_multimin_function * f, const gsl_vector * x, const gsl_vector * step_size)

   The function :func:`gsl_multimin_fdfminimizer_set` initializes the minimizer :data:`s` to minimize the function
   :data:`fdf` starting from the initial point :data:`x`.  The size of the
   first trial step is given by :data:`step_size`.  The accuracy of the line
   minimization is specified by :data:`tol`.  The precise meaning of this
   parameter depends on the method used.  Typically the line minimization
   is considered successful if the gradient of the function :math:`g` is
   orthogonal to the current search direction :math:`p` to a relative
   accuracy of :data:`tol`, where :math:`p \cdot g < tol |p| |g|`.
   A :data:`tol` value of 0.1 is 
   suitable for most purposes, since line minimization only needs to
   be carried out approximately.    Note that setting :data:`tol` to zero will
   force the use of "exact" line-searches, which are extremely expensive.

   The function :func:`gsl_multimin_fminimizer_set` initializes the minimizer :data:`s` to minimize the function
   :data:`f`, starting from the initial point
   :data:`x`. The size of the initial trial steps is given in vector
   :data:`step_size`. The precise meaning of this parameter depends on the
   method used. 

.. function:: void gsl_multimin_fdfminimizer_free (gsl_multimin_fdfminimizer * s)
              void gsl_multimin_fminimizer_free (gsl_multimin_fminimizer * s)

   This function frees all the memory associated with the minimizer
   :data:`s`.

.. function:: const char * gsl_multimin_fdfminimizer_name (const gsl_multimin_fdfminimizer * s)
              const char * gsl_multimin_fminimizer_name (const gsl_multimin_fminimizer * s)

   This function returns a pointer to the name of the minimizer.  For example::

      printf ("s is a '%s' minimizer\n", gsl_multimin_fdfminimizer_name (s));

   would print something like :code:`s is a 'conjugate_pr' minimizer`.

Providing a function to minimize
================================

You must provide a parametric function of :math:`n` variables for the
minimizers to operate on.  You may also need to provide a routine which
calculates the gradient of the function and a third routine which
calculates both the function value and the gradient together.  In order
to allow for general parameters the functions are defined by the
following data types:

.. type:: gsl_multimin_function_fdf

   This data type defines a general function of :math:`n` variables with
   parameters and the corresponding gradient vector of derivatives,

   :code:`double (* f) (const gsl_vector * x, void * params)`

      this function should return the result
      :math:`f(x,params)` for argument :data:`x` and parameters :data:`params`.
      If the function cannot be computed, an error value of :macro:`GSL_NAN`
      should be returned.

   :code:`void (* df) (const gsl_vector * x, void * params, gsl_vector * g)`

      this function should store the :data:`n`-dimensional gradient

      .. only:: not texinfo

         .. math:: g_i = \partial f(x,\hbox{\it params}) / \partial x_i

      .. only:: texinfo

         ::

            g_i = d f(x,params) / d x_i
            
      in the vector :data:`g` for argument :data:`x` 
      and parameters :data:`params`, returning an appropriate error code if the
      function cannot be computed.

   :code:`void (* fdf) (const gsl_vector * x, void * params, double * f, gsl_vector * g)`

      This function should set the values of the :data:`f` and :data:`g` as above,
      for arguments :data:`x` and parameters :data:`params`.  This function
      provides an optimization of the separate functions for :math:`f(x)` and
      :math:`g(x)`---it is always faster to compute the function and its
      derivative at the same time.

   :code:`size_t n`

      the dimension of the system, i.e. the number of components of the
      vectors :data:`x`.

   :code:`void * params`

      a pointer to the parameters of the function.

.. type:: gsl_multimin_function

   This data type defines a general function of :math:`n` variables with
   parameters,

   :code:`double (* f) (const gsl_vector * x, void * params)`

      this function should return the result
      :math:`f(x,params)` for argument :data:`x` and parameters :data:`params`.
      If the function cannot be computed, an error value of :macro:`GSL_NAN`
      should be returned.

   :code:`size_t n`

      the dimension of the system, i.e. the number of components of the
      vectors :data:`x`.

   :code:`void * params`

      a pointer to the parameters of the function.

.. _multimin-paraboloid:

The following example function defines a simple two-dimensional
paraboloid with five parameters,

.. include:: examples/multiminfn.c
   :code:

The function can be initialized using the following code::

  gsl_multimin_function_fdf my_func;

  /* Paraboloid center at (1,2), scale factors (10, 20), 
     minimum value 30 */
  double p[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 }; 

  my_func.n = 2;  /* number of function components */
  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.params = (void *)p;

Iteration
=========

The following function drives the iteration of each algorithm.  The
function performs one iteration to update the state of the minimizer.
The same function works for all minimizers so that different methods can
be substituted at runtime without modifications to the code.

.. function:: int gsl_multimin_fdfminimizer_iterate (gsl_multimin_fdfminimizer * s)
              int gsl_multimin_fminimizer_iterate (gsl_multimin_fminimizer * s)

   These functions perform a single iteration of the minimizer :data:`s`.
   If the iteration encounters an unexpected problem then an error code
   will be returned.  The error code :macro:`GSL_ENOPROG` signifies that
   the minimizer is unable to improve on its current estimate, either due
   to numerical difficulty or because a genuine local minimum has been
   reached.

The minimizer maintains a current best estimate of the minimum at all
times.  This information can be accessed with the following auxiliary
functions,

.. function:: gsl_vector * gsl_multimin_fdfminimizer_x (const gsl_multimin_fdfminimizer * s)
              gsl_vector * gsl_multimin_fminimizer_x (const gsl_multimin_fminimizer * s)
              double gsl_multimin_fdfminimizer_minimum (const gsl_multimin_fdfminimizer * s)
              double gsl_multimin_fminimizer_minimum (const gsl_multimin_fminimizer * s)
              gsl_vector * gsl_multimin_fdfminimizer_gradient (const gsl_multimin_fdfminimizer * s)
              gsl_vector * gsl_multimin_fdfminimizer_dx (const gsl_multimin_fdfminimizer * s)
              double gsl_multimin_fminimizer_size (const gsl_multimin_fminimizer * s)

   These functions return the current best estimate of the location of the
   minimum, the value of the function at that point, its gradient, the last
   step increment of the estimate, and minimizer specific characteristic size for the minimizer :data:`s`.

.. function:: int gsl_multimin_fdfminimizer_restart (gsl_multimin_fdfminimizer * s)

   This function resets the minimizer :data:`s` to use the current point as a
   new starting point.

Stopping Criteria
=================

A minimization procedure should stop when one of the following
conditions is true:

* A minimum has been found to within the user-specified precision.
* A user-specified maximum number of iterations has been reached.
* An error has occurred.

The handling of these conditions is under user control.  The functions
below allow the user to test the precision of the current result.

.. function:: int gsl_multimin_test_gradient (const gsl_vector * g, double epsabs)

   This function tests the norm of the gradient :data:`g` against the
   absolute tolerance :data:`epsabs`. The gradient of a multidimensional
   function goes to zero at a minimum. The test returns :macro:`GSL_SUCCESS`
   if the following condition is achieved,

   .. only:: not texinfo

      .. math:: |g| < \hbox{\it epsabs}

   .. only:: texinfo

      ::

         |g| < epsabs

   and returns :macro:`GSL_CONTINUE` otherwise.  A suitable choice of
   :data:`epsabs` can be made from the desired accuracy in the function for
   small variations in :math:`x`.  The relationship between these quantities
   is given by :math:`\delta{f} = g\,\delta{x}`.

.. function:: int gsl_multimin_test_size (const double size, double epsabs)

   This function tests the minimizer specific characteristic
   size (if applicable to the used minimizer) against absolute tolerance :data:`epsabs`. 
   The test returns :macro:`GSL_SUCCESS` if the size is smaller than tolerance,
   otherwise :macro:`GSL_CONTINUE` is returned.

Algorithms with Derivatives
===========================

There are several minimization methods available. The best choice of
algorithm depends on the problem.  The algorithms described in this
section use the value of the function and its gradient at each
evaluation point.

.. type:: gsl_multimin_fdfminimizer_type

   This type specifies a minimization algorithm using gradients.

   .. index::
      single: Fletcher-Reeves conjugate gradient algorithm, minimization
      single: Conjugate gradient algorithm, minimization
      single: minimization, conjugate gradient algorithm

   .. var:: gsl_multimin_fdfminimizer_type * gsl_multimin_fdfminimizer_conjugate_fr

      This is the Fletcher-Reeves conjugate gradient algorithm. The conjugate
      gradient algorithm proceeds as a succession of line minimizations. The
      sequence of search directions is used to build up an approximation to the
      curvature of the function in the neighborhood of the minimum.  

      An initial search direction :data:`p` is chosen using the gradient, and line
      minimization is carried out in that direction.  The accuracy of the line
      minimization is specified by the parameter :data:`tol`.  The minimum
      along this line occurs when the function gradient :data:`g` and the search direction
      :data:`p` are orthogonal.  The line minimization terminates when
      :math:`p\cdot g < tol |p| |g|`. The
      search direction is updated  using the Fletcher-Reeves formula
      :math:`p' = g' - \beta p` where :math:`\beta=-|g'|^2/|g|^2`, and
      the line minimization is then repeated for the new search
      direction.

   .. index::
      single: Polak-Ribiere algorithm, minimization
      single: minimization, Polak-Ribiere algorithm

   .. var:: gsl_multimin_fdfminimizer_type * gsl_multimin_fdfminimizer_conjugate_pr

      This is the Polak-Ribiere conjugate gradient algorithm.  It is similar
      to the Fletcher-Reeves method, differing only in the choice of the
      coefficient :math:`\beta`. Both methods work well when the evaluation
      point is close enough to the minimum of the objective function that it
      is well approximated by a quadratic hypersurface.

   .. index::
      single: BFGS algorithm, minimization
      single: minimization, BFGS algorithm

   .. var:: gsl_multimin_fdfminimizer_type * gsl_multimin_fdfminimizer_vector_bfgs2
            gsl_multimin_fdfminimizer_type * gsl_multimin_fdfminimizer_vector_bfgs

      These methods use the vector Broyden-Fletcher-Goldfarb-Shanno (BFGS)
      algorithm.  This is a quasi-Newton method which builds up an approximation
      to the second derivatives of the function :math:`f` using the difference
      between successive gradient vectors.  By combining the first and second
      derivatives the algorithm is able to take Newton-type steps towards the
      function minimum, assuming quadratic behavior in that region.

      The :code:`bfgs2` version of this minimizer is the most efficient
      version available, and is a faithful implementation of the line
      minimization scheme described in Fletcher's *Practical Methods of
      Optimization*, Algorithms 2.6.2 and 2.6.4.  It supersedes the original
      :code:`bfgs` routine and requires substantially fewer function and
      gradient evaluations.  The user-supplied tolerance :data:`tol`
      corresponds to the parameter :math:`\sigma` used by Fletcher.  A value
      of 0.1 is recommended for typical use (larger values correspond to
      less accurate line searches).

   .. index::
      single: steepest descent algorithm, minimization
      single: minimization, steepest descent algorithm

   .. var:: gsl_multimin_fdfminimizer_type * gsl_multimin_fdfminimizer_steepest_descent

      The steepest descent algorithm follows the downhill gradient of the
      function at each step. When a downhill step is successful the step-size
      is increased by a factor of two.  If the downhill step leads to a higher
      function value then the algorithm backtracks and the step size is
      decreased using the parameter :data:`tol`.  A suitable value of :data:`tol`
      for most applications is 0.1.  The steepest descent method is
      inefficient and is included only for demonstration purposes.

Algorithms without Derivatives
==============================

The algorithms described in this section use only the value of the function
at each evaluation point.

.. type:: gsl_multimin_fminimizer_type

   This type specifies minimization algorithms which do not use gradients.

   .. index::
      single: Nelder-Mead simplex algorithm for minimization
      single: simplex algorithm, minimization
      single: minimization, simplex algorithm

   .. var:: gsl_multimin_fminimizer_type * gsl_multimin_fminimizer_nmsimplex2
            gsl_multimin_fminimizer_type * gsl_multimin_fminimizer_nmsimplex

      These methods use the Simplex algorithm of Nelder and Mead. 
      Starting from the initial vector :math:`x = p_0`, the algorithm
      constructs an additional :math:`n` vectors :math:`p_i`
      using the step size vector :math:`s = step\_size`
      as follows:

      .. only:: not texinfo

         .. math::

            p_0 & = (x_0, x_1, \cdots , x_n) \\
            p_1 & = (x_0 + s_0, x_1, \cdots , x_n) \\
            p_2 & = (x_0, x_1 + s_1, \cdots , x_n) \\
            \dots &= \dots \\
            p_n & = (x_0, x_1, \cdots , x_n + s_n)

      .. only:: texinfo

         ::

            p_0 = (x_0, x_1, ... , x_n) 
            p_1 = (x_0 + s_0, x_1, ... , x_n) 
            p_2 = (x_0, x_1 + s_1, ... , x_n) 
            ... = ...
            p_n = (x_0, x_1, ... , x_n + s_n)

      These vectors form the :math:`n+1` vertices of a simplex in :math:`n`
      dimensions.  On each iteration the algorithm uses simple geometrical
      transformations to update the vector corresponding to the highest
      function value.  The geometric transformations are reflection,
      reflection followed by expansion, contraction and multiple
      contraction.  Using these transformations the simplex moves through
      the space towards the minimum, where it contracts itself.

      After each iteration, the best vertex is returned.  Note, that due to
      the nature of the algorithm not every step improves the current
      best parameter vector.  Usually several iterations are required.

      The minimizer-specific characteristic size is calculated as the
      average distance from the geometrical center of the simplex to all its
      vertices.  This size can be used as a stopping criteria, as the
      simplex contracts itself near the minimum. The size is returned by the
      function :func:`gsl_multimin_fminimizer_size`.

      The :type:`gsl_multimin_fminimizer_nmsimplex2` version of this minimiser is
      a new :math:`O(N)` operations
      implementation of the earlier :math:`O(N^2)` operations
      :type:`gsl_multimin_fminimizer_nmsimplex`
      minimiser.  It uses the same underlying algorithm, but the simplex
      updates are computed more efficiently for high-dimensional problems.
      In addition, the size of simplex is calculated as the RMS
      distance of each vertex from the center rather than the mean distance,
      allowing a linear update of this quantity on each step.  The memory usage is
      :math:`O(N^2)` for both algorithms.

   .. var:: gsl_multimin_fminimizer_type * gsl_multimin_fminimizer_nmsimplex2rand

      This method is a variant of :type:`gsl_multimin_fminimizer_nmsimplex2` which initialises the
      simplex around the starting point :data:`x` using a randomly-oriented
      set of basis vectors instead of the fixed coordinate axes. The
      final dimensions of the simplex are scaled along the coordinate axes by the
      vector :data:`step_size`.  The randomisation uses a simple deterministic
      generator so that repeated calls to :func:`gsl_multimin_fminimizer_set` for
      a given solver object will vary the orientation in a well-defined way.

Examples
========

This example program finds the minimum of the :ref:`paraboloid function <multimin-paraboloid>`
defined earlier.  The location of the minimum is offset from the origin
in :math:`x` and :math:`y`, and the function value at the minimum is
non-zero. The main program is given below, it requires the example
function given earlier in this chapter.

.. include:: examples/multimin.c
   :code:

The initial step-size is chosen as 0.01, a conservative estimate in this
case, and the line minimization parameter is set at 0.0001.  The program
terminates when the norm of the gradient has been reduced below
0.001. The output of the program is shown below,

.. include:: examples/multimin.txt
   :code:

Note that the algorithm gradually increases the step size as it
successfully moves downhill, as can be seen by plotting the successive
points in :numref:`fig-multimin`.

.. _fig-multimin:

.. figure:: /images/multimin.png
   :scale: 60%

   Function contours with path taken by minimization algorithm

The conjugate gradient algorithm finds the minimum on its second
direction because the function is purely quadratic. Additional
iterations would be needed for a more complicated function.

Here is another example using the Nelder-Mead Simplex algorithm to
minimize the same example object function, as above.

.. include:: examples/nmsimplex.c
   :code:

The minimum search stops when the Simplex size drops to 0.01. The output is
shown below.

.. include:: examples/nmsimplex.txt
   :code:

The simplex size first increases, while the simplex moves towards the
minimum. After a while the size begins to decrease as the simplex
contracts around the minimum.

References and Further Reading
==============================

The conjugate gradient and BFGS methods are described in detail in the
following book,

* R. Fletcher,
  *Practical Methods of Optimization (Second Edition)* Wiley
  (1987), ISBN 0471915475.

A brief description of multidimensional minimization algorithms and
more recent references can be found in,

* C.W. Ueberhuber,
  *Numerical Computation (Volume 2)*, Chapter 14, Section 4.4
  "Minimization Methods", p.: 325--335, Springer (1997), ISBN
  3-540-62057-5.

The simplex algorithm is described in the following paper, 

* J.A. Nelder and R. Mead,
  *A simplex method for function minimization*, Computer Journal
  vol.: 7 (1965), 308--313.
