.. index::
   single: simulated annealing
   single: combinatorial optimization
   single: optimization, combinatorial
   single: energy function
   single: cost function

*******************
Simulated Annealing
*******************

Stochastic search techniques are used when the structure of a space is
not well understood or is not smooth, so that techniques like Newton's
method (which requires calculating Jacobian derivative matrices) cannot
be used. In particular, these techniques are frequently used to solve
combinatorial optimization problems, such as the traveling salesman
problem.

The goal is to find a point in the space at which a real valued
*energy function* (or *cost function*) is minimized.  Simulated
annealing is a minimization technique which has given good results in
avoiding local minima; it is based on the idea of taking a random walk
through the space at successively lower temperatures, where the
probability of taking a step is given by a Boltzmann distribution.

The functions described in this chapter are declared in the header file
:file:`gsl_siman.h`.

Simulated Annealing algorithm
=============================

The simulated annealing algorithm takes random walks through the problem
space, looking for points with low energies; in these random walks, the
probability of taking a step is determined by the Boltzmann distribution,

.. math:: p = e^{-(E_{i+1} - E_i)/(kT)}

if 
:math:`E_{i+1} > E_i`, and
:math:`p = 1` when 
:math:`E_{i+1} \le E_i`.

In other words, a step will occur if the new energy is lower.  If
the new energy is higher, the transition can still occur, and its
likelihood is proportional to the temperature :math:`T` and inversely
proportional to the energy difference 
:math:`E_{i+1} - E_i`.

.. index::
   single: cooling schedule
   single: schedule, cooling

The temperature :math:`T` is initially set to a high value, and a random
walk is carried out at that temperature.  Then the temperature is
lowered very slightly according to a *cooling schedule*, for
example: :math:`T \rightarrow T/\mu_T`
where :math:`\mu_T` is slightly greater than 1. 

The slight probability of taking a step that gives higher energy is what
allows simulated annealing to frequently get out of local minima.

Simulated Annealing functions
=============================

.. function:: void gsl_siman_solve (const gsl_rng * r, void * x0_p, gsl_siman_Efunc_t Ef, gsl_siman_step_t take_step, gsl_siman_metric_t distance, gsl_siman_print_t print_position, gsl_siman_copy_t copyfunc, gsl_siman_copy_construct_t copy_constructor, gsl_siman_destroy_t destructor, size_t element_size, gsl_siman_params_t params)

   This function performs a simulated annealing search through a given
   space.  The space is specified by providing the functions :data:`Ef` and
   :data:`distance`.  The simulated annealing steps are generated using the
   random number generator :data:`r` and the function :data:`take_step`.

   The starting configuration of the system should be given by :data:`x0_p`.
   The routine offers two modes for updating configurations, a fixed-size
   mode and a variable-size mode.  In the fixed-size mode the configuration
   is stored as a single block of memory of size :data:`element_size`.
   Copies of this configuration are created, copied and destroyed
   internally using the standard library functions :func:`malloc`,
   :func:`memcpy` and :func:`free`.  The function pointers :data:`copyfunc`,
   :data:`copy_constructor` and :data:`destructor` should be null pointers in
   fixed-size mode.  In the variable-size mode the functions
   :data:`copyfunc`, :data:`copy_constructor` and :data:`destructor` are used to
   create, copy and destroy configurations internally.  The variable
   :data:`element_size` should be zero in the variable-size mode.

   The :data:`params` structure (described below) controls the run by
   providing the temperature schedule and other tunable parameters to the
   algorithm.

   On exit the best result achieved during the search is placed in
   :data:`x0_p`.  If the annealing process has been successful this
   should be a good approximation to the optimal point in the space.

   If the function pointer :data:`print_position` is not null, a debugging
   log will be printed to :code:`stdout` with the following columns::

      #-iter  #-evals  temperature  position  energy  best_energy

   and the output of the function :data:`print_position` itself.  If
   :data:`print_position` is null then no information is printed.

The simulated annealing routines require several user-specified
functions to define the configuration space and energy function.  The
prototypes for these functions are given below.

.. type:: gsl_siman_Efunc_t

   This function type should return the energy of a configuration :data:`xp`::

      double (*gsl_siman_Efunc_t) (void *xp)

.. type:: gsl_siman_step_t

   This function type should modify the configuration :data:`xp` using a random step
   taken from the generator :data:`r`, up to a maximum distance of
   :data:`step_size`::

      void (*gsl_siman_step_t) (const gsl_rng *r, void *xp, 
                                double step_size)

.. type:: gsl_siman_metric_t

   This function type should return the distance between two configurations
   :data:`xp` and :data:`yp`::

      double (*gsl_siman_metric_t) (void *xp, void *yp)

.. type:: gsl_siman_print_t

   This function type should print the contents of the configuration :data:`xp`::

      void (*gsl_siman_print_t) (void *xp)

.. type:: gsl_siman_copy_t

   This function type should copy the configuration :data:`source` into :data:`dest`::

      void (*gsl_siman_copy_t) (void *source, void *dest)

.. type:: gsl_siman_copy_construct_t

   This function type should create a new copy of the configuration :data:`xp`::

      void * (*gsl_siman_copy_construct_t) (void *xp)

.. type:: gsl_siman_destroy_t

   This function type should destroy the configuration :data:`xp`, freeing its
   memory::

      void (*gsl_siman_destroy_t) (void *xp)

.. type:: gsl_siman_params_t

   These are the parameters that control a run of :func:`gsl_siman_solve`.
   This structure contains all the information needed to control the
   search, beyond the energy function, the step function and the initial
   guess.

   ========================================= ============================================================
   :code:`int n_tries`                       The number of points to try for each step.
   :code:`int iters_fixed_T`                 The number of iterations at each temperature.
   :code:`double step_size`                  The maximum step size in the random walk.
   :code:`double k, t_initial, mu_t, t_min`  The parameters of the Boltzmann distribution and cooling
                                             schedule.
   ========================================= ============================================================

Examples
========

The simulated annealing package is clumsy, and it has to be because it
is written in C, for C callers, and tries to be polymorphic at the same
time.  But here we provide some examples which can be pasted into your
application with little change and should make things easier.

Trivial example
---------------

The first example, in one dimensional Cartesian space, sets up an energy
function which is a damped sine wave; this has many local minima, but
only one global minimum, somewhere between 1.0 and 1.5.  The initial
guess given is 15.5, which is several local minima away from the global
minimum.

.. include:: examples/siman.c
   :code:

:numref:`fig_siman-test` is generated by running
:code:`siman_test` in the following way::

  $ ./siman_test | awk '!/^#/ {print $1, $4}' 
   | graph -y 1.34 1.4 -W0 -X generation -Y position 
   | plot -Tps > siman-test.eps

:numref:`fig_siman-energy` is generated by running
:code:`siman_test` in the following way::

  $ ./siman_test | awk '!/^#/ {print $1, $5}' 
   | graph -y -0.88 -0.83 -W0 -X generation -Y energy 
   | plot -Tps > siman-energy.eps

.. _fig_siman-test:

.. figure:: /images/siman-test.png
   :scale: 60%

   Example of a simulated annealing run: at higher temperatures (early in
   the plot) you see that the solution can fluctuate, but at lower
   temperatures it converges.

.. _fig_siman-energy:

.. figure:: /images/siman-energy.png
   :scale: 60%

   Simulated annealing energy vs generation

.. index::
   single: TSP
   single: traveling salesman problem

Traveling Salesman Problem
--------------------------

The TSP (*Traveling Salesman Problem*) is the classic combinatorial
optimization problem.  I have provided a very simple version of it,
based on the coordinates of twelve cities in the southwestern United
States.  This should maybe be called the *Flying Salesman Problem*,
since I am using the great-circle distance between cities, rather than
the driving distance.  Also: I assume the earth is a sphere, so I don't
use geoid distances.

The :func:`gsl_siman_solve` routine finds a route which is 3490.62
Kilometers long; this is confirmed by an exhaustive search of all
possible routes with the same initial city.

The full code is given below.

.. include:: examples/siman_tsp.c
   :code:

Below are some plots generated in the following way::

  $ ./siman_tsp > tsp.output
  $ grep -v "^#" tsp.output  
   | awk '{print $1, $NF}'
   | graph -y 3300 6500 -W0 -X generation -Y distance 
      -L "TSP - 12 southwest cities"
   | plot -Tps > 12-cities.eps
  $ grep initial_city_coord tsp.output 
    | awk '{print $2, $3}' 
    | graph -X "longitude (- means west)" -Y "latitude" 
       -L "TSP - initial-order" -f 0.03 -S 1 0.1 
    | plot -Tps > initial-route.eps
  $ grep final_city_coord tsp.output 
    | awk '{print $2, $3}' 
    | graph -X "longitude (- means west)" -Y "latitude" 
       -L "TSP - final-order" -f 0.03 -S 1 0.1 
    | plot -Tps > final-route.eps

This is the output showing the initial order of the cities; longitude is
negative, since it is west and I want the plot to look like a map::

  # initial coordinates of cities (longitude and latitude)
  ###initial_city_coord: -105.95 35.68 Santa Fe
  ###initial_city_coord: -112.07 33.54 Phoenix
  ###initial_city_coord: -106.62 35.12 Albuquerque
  ###initial_city_coord: -103.2 34.41 Clovis
  ###initial_city_coord: -107.87 37.29 Durango
  ###initial_city_coord: -96.77 32.79 Dallas
  ###initial_city_coord: -105.92 35.77 Tesuque
  ###initial_city_coord: -107.84 35.15 Grants
  ###initial_city_coord: -106.28 35.89 Los Alamos
  ###initial_city_coord: -106.76 32.34 Las Cruces
  ###initial_city_coord: -108.58 37.35 Cortez
  ###initial_city_coord: -108.74 35.52 Gallup
  ###initial_city_coord: -105.95 35.68 Santa Fe

The optimal route turns out to be::

  # final coordinates of cities (longitude and latitude)
  ###final_city_coord: -105.95 35.68 Santa Fe
  ###final_city_coord: -103.2 34.41 Clovis
  ###final_city_coord: -96.77 32.79 Dallas
  ###final_city_coord: -106.76 32.34 Las Cruces
  ###final_city_coord: -112.07 33.54 Phoenix
  ###final_city_coord: -108.74 35.52 Gallup
  ###final_city_coord: -108.58 37.35 Cortez
  ###final_city_coord: -107.87 37.29 Durango
  ###final_city_coord: -107.84 35.15 Grants
  ###final_city_coord: -106.62 35.12 Albuquerque
  ###final_city_coord: -106.28 35.89 Los Alamos
  ###final_city_coord: -105.92 35.77 Tesuque
  ###final_city_coord: -105.95 35.68 Santa Fe

.. figure:: /images/siman-initial-route.png
   :scale: 60%

   Initial route for the 12 southwestern cities Flying Salesman Problem.

.. figure:: /images/siman-final-route.png
   :scale: 60%

   Final (optimal) route for the 12 southwestern cities Flying Salesman Problem.

Here's a plot of the cost function (energy) versus generation (point in
the calculation at which a new temperature is set) for this problem:

.. figure:: /images/siman-12-cities.png
   :scale: 60%

   Example of a simulated annealing run for the 12 southwestern cities
   Flying Salesman Problem.

References and Further Reading
==============================

Further information is available in the following book,

* *Modern Heuristic Techniques for Combinatorial Problems*, Colin R. Reeves
  (ed.), McGraw-Hill, 1995 (ISBN 0-07-709239-2).
