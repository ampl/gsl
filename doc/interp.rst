.. index:: interpolation, spline

.. _sec_interpolation:

*************
Interpolation
*************

.. include:: include.rst

This chapter describes functions for performing interpolation.  The
library provides a variety of interpolation methods, including Cubic,
Akima, and Steffen splines.  The interpolation types are interchangeable,
allowing different methods to be used without recompiling.
Interpolations can be defined for both normal and periodic boundary
conditions.  Additional functions are available for computing
derivatives and integrals of interpolating functions. Routines
are provided for interpolating both one and two dimensional datasets.

These interpolation methods produce curves that pass through each
datapoint.  To interpolate noisy data with a smoothing curve see
:ref:`chap_basis-splines`.

The functions described in this section are declared in the header files
:file:`gsl_interp.h` and :file:`gsl_spline.h`.

Introduction to 1D Interpolation
================================

Given a set of data points :math:`(x_1, y_1) \dots (x_n, y_n)` the
routines described in this section compute a continuous interpolating
function :math:`y(x)` such that :math:`y(x_i) = y_i`.  The interpolation
is piecewise smooth, and its behavior at the end-points is determined by
the type of interpolation used.

1D Interpolation Functions
==========================

The interpolation function for a given dataset is stored in a
:type:`gsl_interp` object.  These are created by the following functions.

.. type:: gsl_interp

   Workspace for 1D interpolation

.. function:: gsl_interp * gsl_interp_alloc (const gsl_interp_type * T, size_t size)

   This function returns a pointer to a newly allocated interpolation
   object of type :data:`T` for :data:`size` data-points.

.. function:: int gsl_interp_init (gsl_interp * interp, const double xa[], const double ya[], size_t size)

   This function initializes the interpolation object :data:`interp` for the
   data (:data:`xa`, :data:`ya`) where :data:`xa` and :data:`ya` are arrays of size
   :data:`size`.  The interpolation object (:type:`gsl_interp`) does not save
   the data arrays :data:`xa` and :data:`ya` and only stores the static state
   computed from the data.  The :data:`xa` data array is always assumed to be
   strictly ordered, with increasing :math:`x` values; 
   the behavior for other arrangements is not defined.

.. function:: void gsl_interp_free (gsl_interp * interp)

   This function frees the interpolation object :data:`interp`.

1D Interpolation Types
======================

The interpolation library provides the following interpolation types:

.. type:: gsl_interp_type

   .. index:: linear interpolation

   .. var:: gsl_interp_type * gsl_interp_linear

      Linear interpolation.  This interpolation method does not require any
      additional memory.

   .. index:: polynomial interpolation

   .. var:: gsl_interp_type * gsl_interp_polynomial

      Polynomial interpolation.  This method should only be used for
      interpolating small numbers of points because polynomial interpolation
      introduces large oscillations, even for well-behaved datasets.  The
      number of terms in the interpolating polynomial is equal to the number
      of points.

   .. index:: cubic splines

   .. var:: gsl_interp_type * gsl_interp_cspline

      Cubic spline with natural boundary conditions.  The resulting curve is
      piecewise cubic on each interval, with matching first and second
      derivatives at the supplied data-points.  The second derivative is
      chosen to be zero at the first point and last point.

   .. var:: gsl_interp_type * gsl_interp_cspline_periodic

      Cubic spline with periodic boundary conditions.  The resulting curve
      is piecewise cubic on each interval, with matching first and second
      derivatives at the supplied data-points.  The derivatives at the first
      and last points are also matched.  Note that the last point in the
      data must have the same y-value as the first point, otherwise the
      resulting periodic interpolation will have a discontinuity at the
      boundary.

   .. index:: Akima splines

   .. var:: gsl_interp_type * gsl_interp_akima

      Non-rounded Akima spline with natural boundary conditions.  This method
      uses the non-rounded corner algorithm of Wodicka.

   .. var:: gsl_interp_type * gsl_interp_akima_periodic

      Non-rounded Akima spline with periodic boundary conditions.  This method
      uses the non-rounded corner algorithm of Wodicka.

   .. var:: gsl_interp_type * gsl_interp_steffen

      Steffen's method guarantees the monotonicity of the interpolating function
      between the given data points. Therefore, minima and maxima can only occur
      exactly at the data points, and there can never be spurious oscillations
      between data points. The interpolated function is piecewise cubic
      in each interval. The resulting curve and its first derivative
      are guaranteed to be continuous, but the second derivative may be
      discontinuous.

The following related functions are available:

.. function:: const char * gsl_interp_name (const gsl_interp * interp)

   This function returns the name of the interpolation type used by :data:`interp`.
   For example::

      printf ("interp uses '%s' interpolation.\n", gsl_interp_name (interp));

   would print something like::

      interp uses 'cspline' interpolation.

.. function:: unsigned int gsl_interp_min_size (const gsl_interp * interp)
              unsigned int gsl_interp_type_min_size (const gsl_interp_type * T)

   These functions return the minimum number of points required by the
   interpolation object :data:`interp` or interpolation type :data:`T`.  For
   example, Akima spline interpolation requires a minimum of 5 points.

1D Index Look-up and Acceleration
=================================

The state of searches can be stored in a :type:`gsl_interp_accel` object,
which is a kind of iterator for interpolation lookups.

.. type:: gsl_interp_accel

   This workspace stores state variables for interpolation lookups.
   It caches the previous value of an index lookup.  When the subsequent interpolation
   point falls in the same interval its index value can be returned
   immediately.

.. function:: size_t gsl_interp_bsearch (const double x_array[], double x, size_t index_lo, size_t index_hi)

   This function returns the index :math:`i` of the array :data:`x_array` such
   that :code:`x_array[i] <= x < x_array[i+1]`.  The index is searched for
   in the range [:data:`index_lo`, :data:`index_hi`]. |inlinefn|

.. function:: gsl_interp_accel * gsl_interp_accel_alloc (void)

   This function returns a pointer to an accelerator object, which is a
   kind of iterator for interpolation lookups.  It tracks the state of
   lookups, thus allowing for application of various acceleration
   strategies. When multiple interpolants are in use, the same accelerator
   object may be used for all datasets with the same domain
   (:data:`x_array`), but different accelerators should be used for data
   defined on different domains.

.. function:: size_t gsl_interp_accel_find (gsl_interp_accel * a, const double x_array[], size_t size, double x)

   This function performs a lookup action on the data array :data:`x_array`
   of size :data:`size`, using the given accelerator :data:`a`.  This is how
   lookups are performed during evaluation of an interpolation.  The
   function returns an index :math:`i` such that :code:`x_array[i] <= x < x_array[i+1]`.
   |inlinefn|

.. function:: int gsl_interp_accel_reset (gsl_interp_accel * acc);

   This function reinitializes the accelerator object :data:`acc`.  It
   should be used when the cached information is no longer
   applicable---for example, when switching to a new dataset.

.. function:: void gsl_interp_accel_free (gsl_interp_accel* acc)

   This function frees the accelerator object :data:`acc`.

1D Evaluation of Interpolating Functions
========================================

.. function::  double gsl_interp_eval (const gsl_interp * interp, const double xa[], const double ya[], double x, gsl_interp_accel * acc)
               int gsl_interp_eval_e (const gsl_interp * interp, const double xa[], const double ya[], double x, gsl_interp_accel * acc, double * y)

   These functions return the interpolated value of :data:`y` for a given
   point :data:`x`, using the interpolation object :data:`interp`, data
   arrays :data:`xa` and :data:`ya` and the accelerator :data:`acc`.  When
   :data:`x` is outside the range of :data:`xa`, the error code
   :macro:`GSL_EDOM` is returned with a value of :macro:`GSL_NAN` for
   :data:`y`.

.. function:: double gsl_interp_eval_deriv (const gsl_interp * interp, const double xa[], const double ya[], double x, gsl_interp_accel * acc)
              int gsl_interp_eval_deriv_e (const gsl_interp * interp, const double xa[], const double ya[], double x, gsl_interp_accel * acc, double * d)

   These functions return the derivative :data:`d` of an interpolated
   function for a given point :data:`x`, using the interpolation object
   :data:`interp`, data arrays :data:`xa` and :data:`ya` and the accelerator
   :data:`acc`. 

.. function:: double gsl_interp_eval_deriv2 (const gsl_interp * interp, const double xa[], const double ya[], double x, gsl_interp_accel * acc)
              int gsl_interp_eval_deriv2_e (const gsl_interp * interp, const double xa[], const double ya[], double x, gsl_interp_accel * acc, double * d2)

   These functions return the second derivative :data:`d2` of an interpolated
   function for a given point :data:`x`, using the interpolation object
   :data:`interp`, data arrays :data:`xa` and :data:`ya` and the accelerator
   :data:`acc`. 

.. function:: double gsl_interp_eval_integ (const gsl_interp * interp, const double xa[], const double ya[], double a, double b, gsl_interp_accel * acc)
              int gsl_interp_eval_integ_e (const gsl_interp * interp, const double xa[], const double ya[], double a, double b, gsl_interp_accel * acc, double * result)

   These functions return the numerical integral :data:`result` of an
   interpolated function over the range [:data:`a`, :data:`b`], using the
   interpolation object :data:`interp`, data arrays :data:`xa` and :data:`ya` and
   the accelerator :data:`acc`.

1D Higher-level Interface
=========================

The functions described in the previous sections required the user to
supply pointers to the :math:`x` and :math:`y` arrays on each call.  The
following functions are equivalent to the corresponding
:type:`gsl_interp` functions but maintain a copy of this data in the
:type:`gsl_spline` object.  This removes the need to pass both :data:`xa`
and :data:`ya` as arguments on each evaluation. These functions are
defined in the header file :file:`gsl_spline.h`.

.. type:: gsl_spline

   This workspace provides a higher level interface for the
   :type:`gsl_interp` object

.. function:: gsl_spline * gsl_spline_alloc (const gsl_interp_type * T, size_t size)

.. function:: int gsl_spline_init (gsl_spline * spline, const double xa[], const double ya[], size_t size)

.. function:: void gsl_spline_free (gsl_spline * spline)

.. function:: const char * gsl_spline_name (const gsl_spline * spline)

.. function:: unsigned int gsl_spline_min_size (const gsl_spline * spline)

.. function:: double gsl_spline_eval (const gsl_spline * spline, double x, gsl_interp_accel * acc)
              int gsl_spline_eval_e (const gsl_spline * spline, double x, gsl_interp_accel * acc, double * y)

.. function:: double gsl_spline_eval_deriv (const gsl_spline * spline, double x, gsl_interp_accel * acc)
              int gsl_spline_eval_deriv_e (const gsl_spline * spline, double x, gsl_interp_accel * acc, double * d)

.. function:: double gsl_spline_eval_deriv2 (const gsl_spline * spline, double x, gsl_interp_accel * acc)
              int gsl_spline_eval_deriv2_e (const gsl_spline * spline, double x, gsl_interp_accel * acc, double * d2)

.. function:: double gsl_spline_eval_integ (const gsl_spline * spline, double a, double b, gsl_interp_accel * acc)
              int gsl_spline_eval_integ_e (const gsl_spline * spline, double a, double b, gsl_interp_accel * acc, double * result)

1D Interpolation Example Programs
=================================

The following program demonstrates the use of the interpolation and
spline functions.  It computes a cubic spline interpolation of the
10-point dataset :math:`(x_i, y_i)` where :math:`x_i = i + \sin(i)/2` and
:math:`y_i = i + \cos(i^2)` for :math:`i = 0 \dots 9`.

.. include:: examples/interp.c
   :code:

The output is designed to be used with the GNU plotutils
:code:`graph` program::

  $ ./a.out > interp.dat
  $ graph -T ps < interp.dat > interp.ps

.. _fig_interp:

.. figure:: /images/interp.png
   :scale: 60%

   Cubic spline interpolation

:numref:`fig_interp` shows a smooth interpolation of the original points.  The
interpolation method can be changed simply by varying the first argument of
:func:`gsl_spline_alloc`.

The next program demonstrates a periodic cubic spline with 4 data
points.  Note that the first and last points must be supplied with 
the same y-value for a periodic spline.

.. include:: examples/interpp.c
   :code:

The output can be plotted with GNU :code:`graph`::

  $ ./a.out > interp.dat
  $ graph -T ps < interp.dat > interp.ps

.. _fig_interpp:

.. figure:: /images/interpp.png
   :scale: 60%

   Periodic cubic spline interpolation

:numref:`fig_interpp` shows a periodic interpolation of the original points. The
slope of the fitted curve is the same at the beginning and end of the
data, and the second derivative is also.

The next program illustrates the difference between the cubic spline,
Akima, and Steffen interpolation types on a difficult dataset.

.. include:: examples/interp_compare.c
   :code:

.. _fig_interp-compare:

.. figure:: /images/interp_compare.png
   :scale: 60%

   Comparison of different 1D interpolation methods

The output is shown in :numref:`fig_interp-compare`.
The cubic method exhibits a local maxima between the 6th and 7th data points
and continues oscillating for the rest of the data. Akima also shows a
local maxima but recovers and follows the data well after the 7th grid point.
Steffen preserves monotonicity in all intervals and does not exhibit oscillations,
at the expense of having a discontinuous second derivative.

Introduction to 2D Interpolation
================================

Given a set of :math:`x` coordinates :math:`x_1,...,x_m` and a set of
:math:`y` coordinates :math:`y_1,...,y_n`, each in increasing order,
plus a set of function values :math:`z_{ij}`
for each grid point :math:`(x_i,y_j)`, the routines described in this
section compute a continuous interpolation function :math:`z(x,y)` such
that :math:`z(x_i,y_j) = z_{ij}`.

2D Interpolation Functions
==========================

The interpolation function for a given dataset is stored in a
:type:`gsl_interp2d` object. These are created by the following functions.

.. type:: gsl_interp2d

   Workspace for 2D interpolation

.. function:: gsl_interp2d * gsl_interp2d_alloc (const gsl_interp2d_type * T, const size_t xsize, const size_t ysize)

   This function returns a pointer to a newly allocated interpolation
   object of type :data:`T` for :data:`xsize` grid points in the :math:`x`
   direction and :data:`ysize` grid points in the :math:`y` direction.

.. function:: int gsl_interp2d_init (gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const size_t xsize, const size_t ysize)

   This function initializes the interpolation object :data:`interp` for the
   data (:data:`xa`, :data:`ya`, :data:`za`) where :data:`xa` and :data:`ya` are arrays of
   the :math:`x` and :math:`y` grid points of size :data:`xsize` and :data:`ysize`
   respectively, and :data:`za` is an array of function values of size
   :data:`xsize` * :data:`ysize`.  The interpolation object (:type:`gsl_interp2d`) does
   not save the data arrays :data:`xa`, :data:`ya`, and :data:`za` and only stores the
   static state computed from the data. The :data:`xa` and :data:`ya` data arrays
   are always assumed to be strictly ordered, with increasing :math:`x,y` values; 
   the behavior for other arrangements is not defined.

.. function:: void gsl_interp2d_free (gsl_interp2d * interp)

   This function frees the interpolation object :data:`interp`.

2D Interpolation Grids
======================

The 2D interpolation routines access the function values :math:`z_{ij}`
with the following ordering:

.. math:: z_{ij} = za[j*xsize + i]

with :math:`i = 0,...,xsize-1` and :math:`j = 0,...,ysize-1`. However,
for ease of use, the following functions are provided to add and retrieve
elements from the function grid without requiring knowledge of the
internal ordering.

.. function:: int gsl_interp2d_set (const gsl_interp2d * interp, double za[], const size_t i, const size_t j, const double z)

   This function sets the value :math:`z_{ij}` for grid point
   (:data:`i`, :data:`j`) of the array :data:`za` to :data:`z`.

.. function:: double gsl_interp2d_get (const gsl_interp2d * interp, const double za[], const size_t i, const size_t j)

   This function returns the value :math:`z_{ij}` for grid point
   (:data:`i`, :data:`j`) stored in the array :data:`za`.

.. function:: size_t gsl_interp2d_idx (const gsl_interp2d * interp, const size_t i, const size_t j)

   This function returns the index corresponding to the grid point
   (:data:`i`, :data:`j`). The index is given by :math:`j*xsize + i`.

2D Interpolation Types
======================

.. type:: gsl_interp2d_type

   The interpolation library provides the following 2D interpolation types:

   .. index:: bilinear interpolation

   .. var:: gsl_interp2d_type * gsl_interp2d_bilinear

      Bilinear interpolation.  This interpolation method does not require any
      additional memory.

   .. index:: bicubic interpolation

   .. var:: gsl_interp2d_type * gsl_interp2d_bicubic

      Bicubic interpolation.

.. function:: const char * gsl_interp2d_name (const gsl_interp2d * interp)

   This function returns the name of the interpolation type used by :data:`interp`.
   For example::

      printf ("interp uses '%s' interpolation.\n", gsl_interp2d_name (interp));

   would print something like::

      interp uses 'bilinear' interpolation.

.. function:: unsigned int gsl_interp2d_min_size (const gsl_interp2d * interp)
              unsigned int gsl_interp2d_type_min_size (const gsl_interp2d_type * T)

   These functions return the minimum number of points required by the
   interpolation object :data:`interp` or interpolation type :data:`T`.  For
   example, bicubic interpolation requires a minimum of 4 points.

2D Evaluation of Interpolating Functions
========================================

.. function:: double gsl_interp2d_eval (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_interp2d_eval_e (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * z)

   These functions return the interpolated value of :data:`z` for a given
   point (:data:`x`, :data:`y`), using the interpolation object :data:`interp`, data
   arrays :data:`xa`, :data:`ya`, and :data:`za` and the accelerators :data:`xacc`
   and :data:`yacc`.  When :data:`x` is outside the range of :data:`xa` or :data:`y`
   is outside the range of :data:`ya`, the error code
   :macro:`GSL_EDOM` is returned.

.. function:: double gsl_interp2d_eval_extrap (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_interp2d_eval_extrap_e (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * z)

   These functions return the interpolated value of :data:`z` for a given
   point (:data:`x`, :data:`y`), using the interpolation object :data:`interp`, data
   arrays :data:`xa`, :data:`ya`, and :data:`za` and the accelerators :data:`xacc`
   and :data:`yacc`. The functions perform no bounds checking, so
   when :data:`x` is outside the range of :data:`xa` or :data:`y`
   is outside the range of :data:`ya`, extrapolation is performed.

.. function:: double gsl_interp2d_eval_deriv_x (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_interp2d_eval_deriv_x_e (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

   These functions return the interpolated value :data:`d`
   :math:`= \partial z / \partial x` for a given point (:data:`x`, :data:`y`),
   using the interpolation object :data:`interp`, data
   arrays :data:`xa`, :data:`ya`, and :data:`za` and the accelerators :data:`xacc`
   and :data:`yacc`.  When :data:`x` is outside the range of :data:`xa` or :data:`y`
   is outside the range of :data:`ya`, the error code
   :macro:`GSL_EDOM` is returned.

.. function:: double gsl_interp2d_eval_deriv_y (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_interp2d_eval_deriv_y_e (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

   These functions return the interpolated value :data:`d`
   :math:`= \partial z / \partial y` for a given point (:data:`x`, :data:`y`),
   using the interpolation object :data:`interp`, data
   arrays :data:`xa`, :data:`ya`, and :data:`za` and the accelerators :data:`xacc`
   and :data:`yacc`.  When :data:`x` is outside the range of :data:`xa` or :data:`y`
   is outside the range of :data:`ya`, the error code
   :macro:`GSL_EDOM` is returned.

.. function:: double gsl_interp2d_eval_deriv_xx (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_interp2d_eval_deriv_xx_e (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

   These functions return the interpolated value :data:`d`
   :math:`= \partial^2 z / \partial x^2` for a given point (:data:`x`, :data:`y`),
   using the interpolation object :data:`interp`, data
   arrays :data:`xa`, :data:`ya`, and :data:`za` and the accelerators :data:`xacc`
   and :data:`yacc`.  When :data:`x` is outside the range of :data:`xa` or :data:`y`
   is outside the range of :data:`ya`, the error code
   :macro:`GSL_EDOM` is returned.

.. function:: double gsl_interp2d_eval_deriv_yy (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_interp2d_eval_deriv_yy_e (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

   These functions return the interpolated value :data:`d`
   :math:`= \partial^2 z / \partial y^2` for a given point (:data:`x`, :data:`y`),
   using the interpolation object :data:`interp`, data
   arrays :data:`xa`, :data:`ya`, and :data:`za` and the accelerators :data:`xacc`
   and :data:`yacc`.  When :data:`x` is outside the range of :data:`xa` or :data:`y`
   is outside the range of :data:`ya`, the error code
   :macro:`GSL_EDOM` is returned.

.. function:: double gsl_interp2d_eval_deriv_xy (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_interp2d_eval_deriv_xy_e (const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

   These functions return the interpolated value :data:`d`
   :math:`= \partial^2 z / \partial x \partial y` for a given point (:data:`x`, :data:`y`),
   using the interpolation object :data:`interp`, data
   arrays :data:`xa`, :data:`ya`, and :data:`za` and the accelerators :data:`xacc`
   and :data:`yacc`.  When :data:`x` is outside the range of :data:`xa` or :data:`y`
   is outside the range of :data:`ya`, the error code
   :macro:`GSL_EDOM` is returned.

2D Higher-level Interface
=========================

The functions described in the previous sections required the user to
supply pointers to the :math:`x`, :math:`y`, and :math:`z` arrays on each call.
The following functions are equivalent to the corresponding
:code:`gsl_interp2d` functions but maintain a copy of this data in the
:type:`gsl_spline2d` object.  This removes the need to pass :data:`xa`,
:data:`ya`, and :data:`za` as arguments on each evaluation. These functions are
defined in the header file :file:`gsl_spline2d.h`.

.. type:: gsl_spline2d

   This workspace provides a higher level interface for the
   :type:`gsl_interp2d` object

.. function:: gsl_spline2d * gsl_spline2d_alloc (const gsl_interp2d_type * T, size_t xsize, size_t ysize)

.. function:: int gsl_spline2d_init (gsl_spline2d * spline, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize)

.. function:: void gsl_spline2d_free (gsl_spline2d * spline)

.. function:: const char * gsl_spline2d_name (const gsl_spline2d * spline)

.. function:: unsigned int gsl_spline2d_min_size (const gsl_spline2d * spline)

.. function:: double gsl_spline2d_eval (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_spline2d_eval_e (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * z)

.. function:: double gsl_spline2d_eval_extrap (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_spline2d_eval_extrap_e (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * z)

.. function:: double gsl_spline2d_eval_deriv_x (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_spline2d_eval_deriv_x_e (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

.. function:: double gsl_spline2d_eval_deriv_y (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_spline2d_eval_deriv_y_e (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

.. function:: double gsl_spline2d_eval_deriv_xx (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_spline2d_eval_deriv_xx_e (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

.. function:: double gsl_spline2d_eval_deriv_yy (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_spline2d_eval_deriv_yy_e (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

.. function:: double gsl_spline2d_eval_deriv_xy (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc)
              int gsl_spline2d_eval_deriv_xy_e (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d)

.. function:: int gsl_spline2d_set (const gsl_spline2d * spline, double za[], const size_t i, const size_t j, const double z)

.. function:: double gsl_spline2d_get (const gsl_spline2d * spline, const double za[], const size_t i, const size_t j)

   This function returns the value :math:`z_{ij}` for grid point
   (:data:`i`, :data:`j`) stored in the array :data:`za`.

2D Interpolation Example programs
=================================

The following example performs bilinear interpolation on the unit
square, using :math:`z` values of :math:`(0,1,0.5,1)` going clockwise
around the square.

.. include:: examples/interp2d.c
   :code:

The results of the interpolation are shown in :numref:`fig_interp2d`,
where the corners are labeled with their fixed :math:`z` values.

.. _fig_interp2d:

.. figure:: /images/interp2d.png
   :scale: 60%

   2D interpolation example

References and Further Reading
==============================

Descriptions of the interpolation algorithms and further references can
be found in the following publications:

* C.W. Ueberhuber,
  *Numerical Computation (Volume 1), Chapter 9 "Interpolation"*,
  Springer (1997), ISBN 3-540-62058-3.

* D.M. Young, R.T. Gregory,
  *A Survey of Numerical Mathematics (Volume 1), Chapter 6.8*,
  Dover (1988), ISBN 0-486-65691-8.

* M. Steffen,
  *A simple method for monotonic interpolation in one dimension*,
  Astron. Astrophys. 239, 443-450, 1990.
