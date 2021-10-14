Accuracy
========

AMPL wrappers of the GSL functions always use the highest accuracy
available. A GSL function may take a ``mode`` argument that allows
the accuracy of the function to be reduced in order to improve
performance. Currently the ``mode`` argument is always bound to
``GSL_PREC_DOUBLE`` denoting double-precision, a relative accuracy
of approximately :math:`2 \times 10^{-16}`.

.. index::
  accuracy
  GSL_PREC_DOUBLE
