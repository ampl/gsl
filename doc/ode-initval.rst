.. index::
   single: differential equations, initial value problems
   single: initial value problems, differential equations
   single: ordinary differential equations, initial value problem
   single: ODEs, initial value problems

*******************************
Ordinary Differential Equations
*******************************

This chapter describes functions for solving ordinary differential
equation (ODE) initial value problems.  The library provides a variety
of low-level methods, such as Runge-Kutta and Bulirsch-Stoer routines,
and higher-level components for adaptive step-size control. The
components can be combined by the user to achieve the desired
solution, with full access to any intermediate steps. A driver object
can be used as a high level wrapper for easy use of low level
functions.

These functions are declared in the header file :file:`gsl_odeiv2.h`.
This is a new interface in version 1.15 and uses the prefix
:code:`gsl_odeiv2` for all functions.  It is recommended over the
previous :code:`gsl_odeiv` implementation defined in :file:`gsl_odeiv.h`
The old interface has been retained under the original name for
backwards compatibility.

Defining the ODE System
=======================

The routines solve the general :math:`n`-dimensional first-order system,

.. only:: not texinfo

   .. math:: {dy_i(t) \over dt} = f_i (t, y_1(t), \dots y_n(t))

.. only:: texinfo

   ::

      dy_i(t)/dt = f_i(t, y_1(t), ..., y_n(t))

for :math:`i = 1, \dots, n`.  The stepping functions rely on the vector
of derivatives :math:`f_i` and the Jacobian matrix,

.. only:: not texinfo

   .. math:: J_{ij} = \partial f_i(t, y(t)) / \partial y_j

.. only:: texinfo

   ::

      J_{ij} = df_i(t,y(t)) / dy_j

A system of equations is defined using the :type:`gsl_odeiv2_system`
datatype.

.. type:: gsl_odeiv2_system

   This data type defines a general ODE system with arbitrary parameters.

   :code:`int (* function) (double t, const double y[], double dydt[], void * params)`

      This function should store the vector elements
      :math:`f_i(t,y,params)` in the array :data:`dydt`,
      for arguments (:data:`t`, :data:`y`) and parameters :data:`params`.

      The function should return :macro:`GSL_SUCCESS` if the calculation was
      completed successfully. Any other return value indicates an error. A
      special return value :macro:`GSL_EBADFUNC` causes :code:`gsl_odeiv2`
      routines to immediately stop and return. If :code:`function` 
      is modified (for example contents of :data:`params`), the user must call an
      appropriate reset function (:func:`gsl_odeiv2_driver_reset`,
      :func:`gsl_odeiv2_evolve_reset` or :func:`gsl_odeiv2_step_reset`)
      before continuing. Use return values
      distinct from standard GSL error codes to distinguish your function as
      the source of the error.

   .. index::
      single: Jacobian matrix, ODEs

   :code:`int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params)`

      This function should store the vector of derivative elements

      .. only:: not texinfo

         .. math:: \partial f_i(t,y,params) / \partial t
      
      .. only:: texinfo
   
         ::
      
            df_i(t,y,params)/dt

      in the array :data:`dfdt` and the Jacobian matrix :math:`J_{ij}`
      in the array :data:`dfdy`, regarded as a row-ordered
      matrix :code:`J(i,j) = dfdy[i * dimension + j]` where :code:`dimension`
      is the dimension of the system. 

      Not all of the stepper algorithms of :code:`gsl_odeiv2` make use of the
      Jacobian matrix, so it may not be necessary to provide this function
      (the :code:`jacobian` element of the struct can be replaced by a null
      pointer for those algorithms).

      The function should return :macro:`GSL_SUCCESS` if the calculation was
      completed successfully. Any other return value indicates an error. A
      special return value :macro:`GSL_EBADFUNC` causes :code:`gsl_odeiv2`
      routines to immediately stop and return. If :code:`jacobian` 
      is modified (for example contents of :data:`params`), the user must call an
      appropriate reset function (:func:`gsl_odeiv2_driver_reset`,
      :func:`gsl_odeiv2_evolve_reset` or :func:`gsl_odeiv2_step_reset`)
      before continuing. Use return values
      distinct from standard GSL error codes to distinguish your function as
      the source of the error.

   :code:`size_t dimension`

      This is the dimension of the system of equations.

   :code:`void * params`

      This is a pointer to the arbitrary parameters of the system.

Stepping Functions
==================

The lowest level components are the *stepping functions* which
advance a solution from time :math:`t` to :math:`t+h` for a fixed
step-size :math:`h` and estimate the resulting local error.

.. type:: gsl_odeiv2_step

   This contains internal parameters for a stepping function.

.. function:: gsl_odeiv2_step * gsl_odeiv2_step_alloc (const gsl_odeiv2_step_type * T, size_t dim)

   This function returns a pointer to a newly allocated instance of a
   stepping function of type :data:`T` for a system of :data:`dim`
   dimensions. Please note that if you use a stepper method that
   requires access to a driver object, it is advisable to use a driver
   allocation method, which automatically allocates a stepper, too.

.. function:: int gsl_odeiv2_step_reset (gsl_odeiv2_step * s)

   This function resets the stepping function :data:`s`.  It should be used
   whenever the next use of :data:`s` will not be a continuation of a
   previous step.

.. function:: void gsl_odeiv2_step_free (gsl_odeiv2_step * s)

   This function frees all the memory associated with the stepping function
   :data:`s`.

.. function:: const char * gsl_odeiv2_step_name (const gsl_odeiv2_step * s)

   This function returns a pointer to the name of the stepping function.
   For example::

      printf ("step method is '%s'\n", gsl_odeiv2_step_name (s));

   would print something like :code:`step method is 'rkf45'`.

.. function:: unsigned int gsl_odeiv2_step_order (const gsl_odeiv2_step * s)

   This function returns the order of the stepping function on the previous
   step. The order can vary if the stepping function itself is adaptive.

.. function:: int gsl_odeiv2_step_set_driver (gsl_odeiv2_step * s, const gsl_odeiv2_driver * d)

   This function sets a pointer of the driver object :data:`d` for stepper
   :data:`s`, to allow the stepper to access control (and evolve) object
   through the driver object. This is a requirement for some steppers, to
   get the desired error level for internal iteration of
   stepper. Allocation of a driver object calls this function
   automatically.

.. function:: int gsl_odeiv2_step_apply (gsl_odeiv2_step * s, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv2_system * sys)

   This function applies the stepping function :data:`s` to the system of
   equations defined by :data:`sys`, using the step-size :data:`h` to advance
   the system from time :data:`t` and state :data:`y` to time :data:`t` + :data:`h`.
   The new state of the system is stored in :data:`y` on output, with an
   estimate of the absolute error in each component stored in :data:`yerr`.
   If the argument :data:`dydt_in` is not null it should point an array
   containing the derivatives for the system at time :data:`t` on input. This
   is optional as the derivatives will be computed internally if they are
   not provided, but allows the reuse of existing derivative information.
   On output the new derivatives of the system at time :data:`t` + :data:`h` will
   be stored in :data:`dydt_out` if it is not null.

   The stepping function returns :macro:`GSL_FAILURE` if it is unable to
   compute the requested step. Also, if the user-supplied functions
   defined in the system :data:`sys` return a status other than
   :macro:`GSL_SUCCESS` the step will be aborted. In that case, the
   elements of :data:`y` will be restored to their pre-step values and the
   error code from the user-supplied function will be returned. Failure
   may be due to a singularity in the system or too large step-size
   :data:`h`. In that case the step should be attempted again with a
   smaller step-size, e.g. :data:`h` / 2.

   If the driver object is not appropriately set via
   :func:`gsl_odeiv2_step_set_driver` for those steppers that need it, the
   stepping function returns :macro:`GSL_EFAULT`. If the user-supplied
   functions defined in the system :data:`sys` returns :macro:`GSL_EBADFUNC`,
   the function returns immediately with the same return code. In this
   case the user must call :func:`gsl_odeiv2_step_reset` before calling
   this function again.

The following algorithms are available. Please note that algorithms
which use step doubling for error estimation apply the more accurate
values from two half steps instead of values from a single step for
the new state :data:`y`.

.. type:: gsl_odeiv2_step_type

   .. index::
      single: RK2, Runge-Kutta method
      single: Runge-Kutta methods, ordinary differential equations

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_rk2

      Explicit embedded Runge-Kutta (2, 3) method.

   .. index::
      single: RK4, Runge-Kutta method

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_rk4

      Explicit 4th order (classical) Runge-Kutta. Error estimation is
      carried out by the step doubling method. For more efficient estimate
      of the error, use the embedded methods described below.

   .. index::
      single: Fehlberg method, differential equations
      single: RKF45, Runge-Kutta-Fehlberg method

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_rkf45

      Explicit embedded Runge-Kutta-Fehlberg (4, 5) method.  This method is
      a good general-purpose integrator.

   .. index::
      single: Runge-Kutta Cash-Karp method
      single: Cash-Karp, Runge-Kutta method

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_rkck 

      Explicit embedded Runge-Kutta Cash-Karp (4, 5) method.

   .. index::
      single: Runge-Kutta Prince-Dormand method
      single: Prince-Dormand, Runge-Kutta method

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_rk8pd  

      Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.

   .. index:: Implicit Euler method

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_rk1imp 

      Implicit Gaussian first order Runge-Kutta. Also known as implicit
      Euler or backward Euler method. Error estimation is carried out by the
      step doubling method. This algorithm requires the Jacobian and 
      access to the driver object via :func:`gsl_odeiv2_step_set_driver`.

   .. index:: Implicit Runge-Kutta method

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_rk2imp  

      Implicit Gaussian second order Runge-Kutta. Also known as implicit
      mid-point rule. Error estimation is carried out by the step doubling
      method. This stepper requires the Jacobian and access to the driver
      object via :func:`gsl_odeiv2_step_set_driver`.

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_rk4imp  

      Implicit Gaussian 4th order Runge-Kutta. Error estimation is carried
      out by the step doubling method. This algorithm requires the Jacobian
      and access to the driver object via :func:`gsl_odeiv2_step_set_driver`.

   .. index::
      single: Bulirsch-Stoer method
      single: Bader and Deuflhard, Bulirsch-Stoer method.
      single: Deuflhard and Bader, Bulirsch-Stoer method.

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_bsimp   

      Implicit Bulirsch-Stoer method of Bader and Deuflhard. The method is
      generally suitable for stiff problems. This stepper requires the
      Jacobian.

   .. index::
      single: Adams method
      single: multistep methods, ODEs
      single: predictor-corrector method, ODEs
      single: Nordsieck form

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_msadams  

      A variable-coefficient linear multistep Adams method in Nordsieck
      form. This stepper uses explicit Adams-Bashforth (predictor) and
      implicit Adams-Moulton (corrector) methods in :math:`P(EC)^m`
      functional iteration mode. Method order varies dynamically between 1
      and 12. This stepper requires the access to the driver object via
      :func:`gsl_odeiv2_step_set_driver`.

   .. index:: BDF method

   .. var:: gsl_odeiv2_step_type * gsl_odeiv2_step_msbdf

      A variable-coefficient linear multistep backward differentiation
      formula (BDF) method in Nordsieck form. This stepper uses the explicit
      BDF formula as predictor and implicit BDF formula as corrector. A
      modified Newton iteration method is used to solve the system of
      non-linear equations. Method order varies dynamically between 1 and
      5. The method is generally suitable for stiff problems. This stepper
      requires the Jacobian and the access to the driver object via
      :func:`gsl_odeiv2_step_set_driver`.

.. index::
   single: Adaptive step-size control, differential equations

Adaptive Step-size Control
==========================

The control function examines the proposed change to the solution
produced by a stepping function and attempts to determine the optimal
step-size for a user-specified level of error.

.. type:: gsl_odeiv2_control

   This is a workspace for controlling step size.

.. type:: gsl_odeiv2_control_type

   This specifies the control type.

.. function:: gsl_odeiv2_control * gsl_odeiv2_control_standard_new (double eps_abs, double eps_rel, double a_y, double a_dydt)

   The standard control object is a four parameter heuristic based on
   absolute and relative errors :data:`eps_abs` and :data:`eps_rel`, and
   scaling factors :data:`a_y` and :data:`a_dydt` for the system state
   :math:`y(t)` and derivatives :math:`y'(t)` respectively.

   The step-size adjustment procedure for this method begins by computing
   the desired error level :math:`D_i` for each component,

   .. only:: not texinfo

      .. math:: D_i = \epsilon_{abs} + \epsilon_{rel} * (a_{y} |y_i| + a_{dydt} h |y\prime_i|)

   .. only:: texinfo

      ::

         D_i = eps_abs + eps_rel * (a_y |y_i| + a_dydt h |y\prime_i|)

   and comparing it with the observed error :math:`E_i = |yerr_i|`.  If the
   observed error :data:`E` exceeds the desired error level :data:`D` by more
   than 10% for any component then the method reduces the step-size by an
   appropriate factor,

   .. math:: h_{new} = h_{old} * S * (E/D)^{-1/q}

   where :math:`q` is the consistency order of the method (e.g. :math:`q=4` for
   4(5) embedded RK), and :math:`S` is a safety factor of 0.9. The ratio
   :math:`E/D` is taken to be the maximum of the ratios
   :math:`E_i/D_i`. 

   If the observed error :math:`E` is less than 50% of the desired error
   level :data:`D` for the maximum ratio :math:`E_i/D_i` then the algorithm
   takes the opportunity to increase the step-size to bring the error in
   line with the desired level,

   .. math:: h_{new} = h_{old} * S * (E/D)^{-1/(q+1)}

   This encompasses all the standard error scaling methods. To avoid
   uncontrolled changes in the stepsize, the overall scaling factor is
   limited to the range :math:`1/5` to 5.

.. function:: gsl_odeiv2_control * gsl_odeiv2_control_y_new (double eps_abs, double eps_rel)

   This function creates a new control object which will keep the local
   error on each step within an absolute error of :data:`eps_abs` and
   relative error of :data:`eps_rel` with respect to the solution :math:`y_i(t)`.
   This is equivalent to the standard control object with :data:`a_y` = 1 and
   :data:`a_dydt` = 0.

.. function:: gsl_odeiv2_control * gsl_odeiv2_control_yp_new (double eps_abs, double eps_rel)

   This function creates a new control object which will keep the local
   error on each step within an absolute error of :data:`eps_abs` and
   relative error of :data:`eps_rel` with respect to the derivatives of the
   solution :math:`y'_i(t)`.  This is equivalent to the standard control
   object with :data:`a_y` = 0 and :data:`a_dydt` = 1.

.. function:: gsl_odeiv2_control * gsl_odeiv2_control_scaled_new (double eps_abs, double eps_rel, double a_y, double a_dydt, const double scale_abs[], size_t dim)

   This function creates a new control object which uses the same algorithm
   as :func:`gsl_odeiv2_control_standard_new` but with an absolute error
   which is scaled for each component by the array :data:`scale_abs`.
   The formula for :math:`D_i` for this control object is,

   .. only:: not texinfo

      .. math:: D_i = \epsilon_{abs} s_i + \epsilon_{rel} * (a_{y} |y_i| + a_{dydt} h |y\prime_i|)

   .. only:: texinfo

      ::

         D_i = eps_abs * s_i + eps_rel * (a_y |y_i| + a_dydt h |y\prime_i|)

   where :math:`s_i` is the :math:`i`-th component of the array :data:`scale_abs`.
   The same error control heuristic is used by the Matlab ODE suite. 

.. function:: gsl_odeiv2_control * gsl_odeiv2_control_alloc (const gsl_odeiv2_control_type * T)

   This function returns a pointer to a newly allocated instance of a
   control function of type :data:`T`.  This function is only needed for
   defining new types of control functions.  For most purposes the standard
   control functions described above should be sufficient. 

.. function:: int gsl_odeiv2_control_init (gsl_odeiv2_control * c, double eps_abs, double eps_rel, double a_y, double a_dydt)

   This function initializes the control function :data:`c` with the
   parameters :data:`eps_abs` (absolute error), :data:`eps_rel` (relative
   error), :data:`a_y` (scaling factor for y) and :data:`a_dydt` (scaling
   factor for derivatives).

.. function:: void gsl_odeiv2_control_free (gsl_odeiv2_control * c)

   This function frees all the memory associated with the control function
   :data:`c`.

.. function:: int gsl_odeiv2_control_hadjust (gsl_odeiv2_control * c, gsl_odeiv2_step * s, const double y[], const double yerr[], const double dydt[], double * h)

   This function adjusts the step-size :data:`h` using the control function
   :data:`c`, and the current values of :data:`y`, :data:`yerr` and :data:`dydt`.
   The stepping function :data:`step` is also needed to determine the order
   of the method.  If the error in the y-values :data:`yerr` is found to be
   too large then the step-size :data:`h` is reduced and the function returns
   :macro:`GSL_ODEIV_HADJ_DEC`.  If the error is sufficiently small then
   :data:`h` may be increased and :macro:`GSL_ODEIV_HADJ_INC` is returned.  The
   function returns :macro:`GSL_ODEIV_HADJ_NIL` if the step-size is
   unchanged.  The goal of the function is to estimate the largest
   step-size which satisfies the user-specified accuracy requirements for
   the current point.

.. function:: const char * gsl_odeiv2_control_name (const gsl_odeiv2_control * c)

   This function returns a pointer to the name of the control function.
   For example::

      printf ("control method is '%s'\n", gsl_odeiv2_control_name (c));

   would print something like :code:`control method is 'standard'`

.. function:: int gsl_odeiv2_control_errlevel (gsl_odeiv2_control * c, const double y, const double dydt, const double h, const size_t ind, double * errlev)

   This function calculates the desired error level of the :data:`ind`-th component
   to :data:`errlev`. It requires the value (:data:`y`) and value of the derivative
   (:data:`dydt`) of the component, and the current step size :data:`h`.

.. function:: int gsl_odeiv2_control_set_driver (gsl_odeiv2_control * c, const gsl_odeiv2_driver * d)

   This function sets a pointer of the driver object :data:`d` for control
   object :data:`c`.

Evolution
=========

The evolution function combines the results of a stepping function and
control function to reliably advance the solution forward one step
using an acceptable step-size.

.. type:: gsl_odeiv2_evolve

   This workspace contains parameters for controlling the evolution function

.. function:: gsl_odeiv2_evolve * gsl_odeiv2_evolve_alloc (size_t dim)

   This function returns a pointer to a newly allocated instance of an
   evolution function for a system of :data:`dim` dimensions.

.. function:: int gsl_odeiv2_evolve_apply (gsl_odeiv2_evolve * e, gsl_odeiv2_control * con, gsl_odeiv2_step * step, const gsl_odeiv2_system * sys, double * t, double t1, double * h, double y[])

   This function advances the system (:data:`e`, :data:`sys`) from time
   :data:`t` and position :data:`y` using the stepping function :data:`step`.
   The new time and position are stored in :data:`t` and :data:`y` on output.

   The initial step-size is taken as :data:`h`. The control function
   :data:`con` is applied to check whether the local error estimated by the
   stepping function :data:`step` using step-size :data:`h` exceeds the
   required error tolerance. If the error is too high, the step is
   retried by calling :data:`step` with a decreased step-size. This process
   is continued until an acceptable step-size is found. An estimate of
   the local error for the step can be obtained from the components of
   the array :code:`e->yerr[]`.

   If the user-supplied functions defined in the system :data:`sys` returns
   :macro:`GSL_EBADFUNC`, the function returns immediately with the same
   return code. In this case the user must call
   :func:`gsl_odeiv2_step_reset` and
   :func:`gsl_odeiv2_evolve_reset` before calling this function again.

   Otherwise, if the user-supplied functions defined in the system
   :data:`sys` or the stepping function :data:`step` return a status other
   than :macro:`GSL_SUCCESS`, the step is retried with a decreased
   step-size. If the step-size decreases below machine precision, a
   status of :macro:`GSL_FAILURE` is returned if the user functions
   returned :macro:`GSL_SUCCESS`. Otherwise the value returned by user
   function is returned. If no acceptable step can be made, :data:`t` and
   :data:`y` will be restored to their pre-step values and :data:`h` contains
   the final attempted step-size.

   If the step is successful the function returns a suggested step-size
   for the next step in :data:`h`. The maximum time :data:`t1` is guaranteed
   not to be exceeded by the time-step. On the final time-step the value
   of :data:`t` will be set to :data:`t1` exactly.

.. function:: int gsl_odeiv2_evolve_apply_fixed_step (gsl_odeiv2_evolve * e, gsl_odeiv2_control * con, gsl_odeiv2_step * step, const gsl_odeiv2_system * sys, double * t, const double h, double y[])

   This function advances the ODE-system (:data:`e`, :data:`sys`, :data:`con`)
   from time :data:`t` and position :data:`y` using the stepping function
   :data:`step` by a specified step size :data:`h`. If the local error
   estimated by the stepping function exceeds the desired error level,
   the step is not taken and the function returns
   :macro:`GSL_FAILURE`. Otherwise the value returned by user function is
   returned.

.. function:: int gsl_odeiv2_evolve_reset (gsl_odeiv2_evolve * e)

   This function resets the evolution function :data:`e`.  It should be used
   whenever the next use of :data:`e` will not be a continuation of a
   previous step.

.. function:: void gsl_odeiv2_evolve_free (gsl_odeiv2_evolve * e)

   This function frees all the memory associated with the evolution function
   :data:`e`.

.. function:: int gsl_odeiv2_evolve_set_driver (gsl_odeiv2_evolve * e, const gsl_odeiv2_driver * d)

   This function sets a pointer of the driver object :data:`d` for evolve
   object :data:`e`.

.. index::
   single: discontinuities, in ODE systems

If a system has discontinuous changes in the derivatives at known
points, it is advisable to evolve the system between each discontinuity
in sequence.  For example, if a step-change in an external driving
force occurs at times :math:`t_a, t_b` and :math:`t_c` then evolution
should be carried out over the ranges :math:`(t_0,t_a)`,
:math:`(t_a,t_b)`, :math:`(t_b,t_c)`, and :math:`(t_c,t_1)` separately
and not directly over the range :math:`(t_0,t_1)`.

Driver
======

The driver object is a high level wrapper that combines the evolution,
control and stepper objects for easy use.

.. function:: gsl_odeiv2_driver * gsl_odeiv2_driver_alloc_y_new (const gsl_odeiv2_system * sys, const gsl_odeiv2_step_type * T, const double hstart, const double epsabs, const double epsrel)
              gsl_odeiv2_driver * gsl_odeiv2_driver_alloc_yp_new (const gsl_odeiv2_system * sys, const gsl_odeiv2_step_type * T, const double hstart, const double epsabs, const double epsrel)
              gsl_odeiv2_driver * gsl_odeiv2_driver_alloc_standard_new (const gsl_odeiv2_system * sys, const gsl_odeiv2_step_type * T, const double hstart, const double epsabs, const double epsrel, const double a_y, const double a_dydt)
              gsl_odeiv2_driver * gsl_odeiv2_driver_alloc_scaled_new (const gsl_odeiv2_system * sys, const gsl_odeiv2_step_type * T, const double hstart, const double epsabs, const double epsrel, const double a_y, const double a_dydt, const double scale_abs[])

   These functions return a pointer to a newly allocated instance of a
   driver object. The functions automatically allocate and initialise the
   evolve, control and stepper objects for ODE system :data:`sys` using
   stepper type :data:`T`. The initial step size is given in
   :data:`hstart`. The rest of the arguments follow the syntax and
   semantics of the control functions with same name
   (:code:`gsl_odeiv2_control_*_new`).

.. function:: int gsl_odeiv2_driver_set_hmin (gsl_odeiv2_driver * d, const double hmin)

   The function sets a minimum for allowed step size :data:`hmin` for
   driver :data:`d`. Default value is 0.

.. function:: int gsl_odeiv2_driver_set_hmax (gsl_odeiv2_driver * d, const double hmax)

   The function sets a maximum for allowed step size :data:`hmax` for
   driver :data:`d`. Default value is :macro:`GSL_DBL_MAX`.

.. function:: int gsl_odeiv2_driver_set_nmax (gsl_odeiv2_driver * d, const unsigned long int nmax)

   The function sets a maximum for allowed number of steps :data:`nmax` for
   driver :data:`d`. Default value of 0 sets no limit for steps.

.. function:: int gsl_odeiv2_driver_apply (gsl_odeiv2_driver * d, double * t, const double t1, double y[])

   This function evolves the driver system :data:`d` from :data:`t` to
   :data:`t1`. Initially vector :data:`y` should contain the values of
   dependent variables at point :data:`t`. If the function is unable to
   complete the calculation, an error code from
   :func:`gsl_odeiv2_evolve_apply` is returned, and :data:`t` and :data:`y`
   contain the values from last successful step. 

   If maximum number of steps is reached, a value of :macro:`GSL_EMAXITER`
   is returned. If the step size drops below minimum value, the function
   returns with :macro:`GSL_ENOPROG`. If the user-supplied functions
   defined in the system :data:`sys` returns :macro:`GSL_EBADFUNC`, the
   function returns immediately with the same return code. In this case
   the user must call :func:`gsl_odeiv2_driver_reset` before calling this
   function again.

.. function:: int gsl_odeiv2_driver_apply_fixed_step (gsl_odeiv2_driver * d, double * t, const double h, const unsigned long int n, double y[])

   This function evolves the driver system :data:`d` from :data:`t` with
   :data:`n` steps of size :data:`h`. If the function is unable to complete
   the calculation, an error code from
   :func:`gsl_odeiv2_evolve_apply_fixed_step` is returned, and :data:`t` and
   :data:`y` contain the values from last successful step.

.. function:: int gsl_odeiv2_driver_reset (gsl_odeiv2_driver * d)

   This function resets the evolution and stepper objects.

.. function:: int gsl_odeiv2_driver_reset_hstart (gsl_odeiv2_driver * d, const double hstart)

   The routine resets the evolution and stepper objects and sets new
   initial step size to :data:`hstart`. This function can be used e.g. to
   change the direction of integration.

.. function:: int gsl_odeiv2_driver_free (gsl_odeiv2_driver * d)

   This function frees the driver object, and the related evolution,
   stepper and control objects.

Examples
========

.. index::
   single: Van der Pol oscillator, example

The following program solves the second-order nonlinear Van der Pol
oscillator equation,

.. math:: u''(t) + \mu u'(t) (u(t)^2 - 1) + u(t) = 0

This can be converted into a first order system suitable for use with
the routines described in this chapter by introducing a separate
variable for the velocity, :math:`v = u'(t)`,

.. only:: not texinfo

   .. math::

      u' &= v \\
      v' &= -u + \mu v (1-u^2)

.. only:: texinfo

   ::

      u' = v
      v' = -u + \mu v (1-u^2)

The program begins by defining functions for these derivatives and
their Jacobian. The main function uses driver level functions to solve
the problem. The program evolves the solution from :math:`(u, v) = (1, 0)`
at :math:`t = 0` to :math:`t = 100`.  The step-size :math:`h` is
automatically adjusted by the controller to maintain an absolute
accuracy of :math:`10^{-6}`
in the function values :math:`(u, v)`.
The loop in the example prints the solution at the points
:math:`t_i = 1, 2, \dots, 100`.

.. include:: examples/ode-initval.c
   :code:

The user can work with the lower level functions directly, as in
the following example. In this case an intermediate result is printed
after each successful step instead of equidistant time points. 

.. include:: examples/ode-initval-low-level.c
   :code:

For functions with multiple parameters, the appropriate information
can be passed in through the :data:`params` argument in
:type:`gsl_odeiv2_system` definition (:data:`mu` in this example) by using
a pointer to a struct.

.. figure:: /images/ode-vdp.png

   Numerical solution of the Van der Pol oscillator equation 
   using Prince-Dormand 8th order Runge-Kutta.

It is also possible to work with a non-adaptive integrator, using only
the stepping function itself,
:func:`gsl_odeiv2_driver_apply_fixed_step` or
:func:`gsl_odeiv2_evolve_apply_fixed_step`. The following program uses
the driver level function, with fourth-order
Runge-Kutta stepping function with a fixed stepsize of
0.001.

.. include:: examples/odefixed.c
   :code:

References and Further Reading
==============================

* Ascher, U.M., Petzold, L.R., *Computer Methods for Ordinary
  Differential and Differential-Algebraic Equations*, SIAM, 
  Philadelphia, 1998.

* Hairer, E., Norsett, S. P., Wanner, G., *Solving Ordinary Differential 
  Equations I: Nonstiff Problems*, Springer, Berlin, 1993.

* Hairer, E., Wanner, G., *Solving Ordinary Differential 
  Equations II: Stiff and Differential-Algebraic Problems*,
  Springer, Berlin, 1996.

Many of the basic Runge-Kutta formulas can be found in the Handbook of
Mathematical Functions,

* Abramowitz & Stegun (eds.), *Handbook of Mathematical Functions*,
  Section 25.5.

The implicit Bulirsch-Stoer algorithm :code:`bsimp` is described in the
following paper,

* G. Bader and P. Deuflhard, "A Semi-Implicit Mid-Point Rule for Stiff
  Systems of Ordinary Differential Equations.", Numer.: Math.: 41, 373--398,
  1983.

The Adams and BDF multistep methods :code:`msadams` and :code:`msbdf`
are based on the following articles,

* G. D. Byrne and A. C. Hindmarsh, "A Polyalgorithm for the
  Numerical Solution of Ordinary Differential Equations.",
  ACM Trans. Math. Software, 1, 71--96, 1975.

* P. N. Brown, G. D. Byrne and A. C. Hindmarsh, "VODE: A
  Variable-coefficient ODE Solver.", SIAM J. Sci. Stat. Comput. 10,
  1038--1051, 1989.

* A. C. Hindmarsh, P. N. Brown, K. E. Grant, S. L. Lee, R. Serban,
  D. E. Shumaker and C. S. Woodward, "SUNDIALS: Suite of
  Nonlinear and Differential/Algebraic Equation Solvers.", ACM
  Trans. Math. Software 31, 363--396, 2005.
