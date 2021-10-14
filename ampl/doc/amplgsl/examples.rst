Examples
========

Finding Minimum of the Gamma Function
-------------------------------------

The following model (`gamma.mod
<https://raw.github.com/vitaut/ampl/master/models/gsl/gamma.mod>`_) finds a
minimum of the :ref:`gamma function <gamma-functions>` :math:`\Gamma(x)` for
:math:`x > 0`:

.. code-block:: none

  # Find minimum of the gamma function for x > 0.

  include gsl.ampl;

  var x >= 1e-5;
  minimize obj: gsl_sf_gamma(x);
  solve;
  print x;

Solving this models gives an optimal solution :math:`x \approx 1.46163` which
can be verified `here <http://oeis.org/A030169>`_.

Hock and Schittkowski Models
----------------------------

The `hs068.mod <http://ampl.github.io/models/gsl/hs068.mod>`_ and
`hs069.mod <http://ampl.github.io/models/gsl/hs069.mod>`_
models demonstrate the use of :func:`gsl_cdf_ugaussian_P`. They are taken
form the Robert Vanderbei's `collection of nonlinear models
<http://orfe.princeton.edu/~rvdb/ampl/nlmodels/>`_ and adapted for AMPLGSL.
These models can be solved with a nonlinear AMPL solver such as MINOS:

.. code-block:: none

  $ ampl hs068.mod
  obj = -0.261841

  MINOS 5.51: optimal solution found.
  29 iterations, objective -0.9204250037
  Nonlin evals: obj = 57, grad = 56, constrs = 57, Jac = 56.
  x [*] :=
  1  0.0678586
  2  3.64621
  3  0.000266135
  4  0.894855
  ;

  obj = -0.920425

  Best known objective value: -0.920425026
