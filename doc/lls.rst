.. index::
   single: fitting
   single: least squares fit
   single: regression, least squares
   single: weighted linear fits
   single: unweighted linear fits

****************************
Linear Least-Squares Fitting
****************************

This chapter describes routines for performing least squares fits to
experimental data using linear combinations of functions.  The data
may be weighted or unweighted, i.e. with known or unknown errors.  For
weighted data the functions compute the best fit parameters and their
associated covariance matrix.  For unweighted data the covariance
matrix is estimated from the scatter of the points, giving a
variance-covariance matrix.

The functions are divided into separate versions for simple one- or
two-parameter regression and multiple-parameter fits.

.. _sec_lls-overview:

Overview
========

Least-squares fits are found by minimizing :math:`\chi^2`
(chi-squared), the weighted sum of squared residuals over :math:`n`
experimental datapoints :math:`(x_i, y_i)` for the model :math:`Y(c,x)`,

.. math:: \chi^2 = \sum_i w_i (y_i - Y(c, x_i))^2

The :math:`p` parameters of the model are :math:`c = \{c_0, c_1, \dots\}`.
The weight factors :math:`w_i` are given by :math:`w_i = 1/\sigma_i^2`
where :math:`\sigma_i` is the experimental error on the data-point
:math:`y_i`.  The errors are assumed to be
Gaussian and uncorrelated. 
For unweighted data the chi-squared sum is computed without any weight factors. 

.. index::
   single: covariance matrix, linear fits

The fitting routines return the best-fit parameters :math:`c` and their
:math:`p \times p` covariance matrix.  The covariance matrix measures the
statistical errors on the best-fit parameters resulting from the 
errors on the data, :math:`\sigma_i`, and is defined
as

.. only:: not texinfo

   .. math:: C_{ab} = \langle \delta c_a \delta c_b \rangle

.. only:: texinfo

   ::

      C_{ab} = <\delta c_a \delta c_b>

where :math:`\langle \, \rangle`
denotes an average over the Gaussian error distributions of the underlying datapoints.

The covariance matrix is calculated by error propagation from the data
errors :math:`\sigma_i`.  The change in a fitted parameter :math:`\delta c_a`
caused by a small change in the data :math:`\delta y_i` is given
by

.. only:: not texinfo

   .. math:: \delta c_a = \sum_i {\partial c_a \over \partial y_i} \delta y_i

.. only:: texinfo

   ::

      \delta c_a = \sum_i (dc_a/dy_i) \delta y_i

allowing the covariance matrix to be written in terms of the errors on the data,

.. only:: not texinfo

   .. math::

      C_{ab} =  \sum_{i,j} {\partial c_a \over \partial y_i}
                           {\partial c_b \over \partial y_j} 
                           \langle \delta y_i \delta y_j \rangle

.. only:: texinfo

   ::

      C_{ab} = \sum_{i,j} (dc_a/dy_i) (dc_b/dy_j) <\delta y_i \delta y_j>

For uncorrelated data the fluctuations of the underlying datapoints satisfy

.. only:: not texinfo

   .. math:: \langle \delta y_i \delta y_j \rangle = \sigma_i^2 \delta_{ij}

.. only:: texinfo

   ::

      <\delta y_i \delta y_j> = \sigma_i^2 \delta_{ij}

giving a corresponding parameter covariance matrix of

.. only:: not texinfo

   .. math:: C_{ab} = \sum_{i} {1 \over w_i} {\partial c_a \over \partial y_i} {\partial c_b \over \partial y_i} 

.. only:: texinfo

   ::

      C_{ab} = \sum_i (1/w_i) (dc_a/dy_i) (dc_b/dy_i) 

.. index::
   single: variance-covariance matrix, linear fits

When computing the covariance matrix for unweighted data, i.e. data with unknown errors, 
the weight factors :math:`w_i` in this sum are replaced by the single estimate
:math:`w = 1/\sigma^2`, where :math:`\sigma^2` is the computed variance of the
residuals about the best-fit model, :math:`\sigma^2 = \sum (y_i - Y(c,x_i))^2 / (n-p)`.  
This is referred to as the *variance-covariance matrix*.

The standard deviations of the best-fit parameters are given by the
square root of the corresponding diagonal elements of
the covariance matrix, :math:`\sigma_{c_a} = \sqrt{C_{aa}}`.
The correlation coefficient of the fit parameters :math:`c_a` and :math:`c_b`
is given by :math:`\rho_{ab} = C_{ab} / \sqrt{C_{aa} C_{bb}}`.

.. index:: linear regression

Linear regression
=================

The functions in this section are used to fit simple one or two
parameter linear regression models. The functions are declared in
the header file :file:`gsl_fit.h`.

Linear regression with a constant term
--------------------------------------
The functions described in this section can be used to perform
least-squares fits to a straight line model, :math:`Y(c,x) = c_0 + c_1 x`.

.. index::
   single: covariance matrix, from linear regression

.. function:: int gsl_fit_linear (const double * x, const size_t xstride, const double * y, const size_t ystride, size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * sumsq)

   This function computes the best-fit linear regression coefficients
   (:data:`c0`, :data:`c1`) of the model :math:`Y = c_0 + c_1 X` for the dataset
   (:data:`x`, :data:`y`), two vectors of length :data:`n` with strides
   :data:`xstride` and :data:`ystride`.  The errors on :data:`y` are assumed unknown so 
   the variance-covariance matrix for the
   parameters (:data:`c0`, :data:`c1`) is estimated from the scatter of the
   points around the best-fit line and returned via the parameters
   (:data:`cov00`, :data:`cov01`, :data:`cov11`).   
   The sum of squares of the residuals from the best-fit line is returned
   in :data:`sumsq`.  Note: the correlation coefficient of the data can be computed using
   :func:`gsl_stats_correlation`, it does not depend on the fit.

.. function:: int gsl_fit_wlinear (const double * x, const size_t xstride, const double * w, const size_t wstride, const double * y, const size_t ystride, size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * chisq)

   This function computes the best-fit linear regression coefficients
   (:data:`c0`, :data:`c1`) of the model :math:`Y = c_0 + c_1 X` for the weighted
   dataset (:data:`x`, :data:`y`), two vectors of length :data:`n` with strides
   :data:`xstride` and :data:`ystride`.  The vector :data:`w`, of length :data:`n`
   and stride :data:`wstride`, specifies the weight of each datapoint. The
   weight is the reciprocal of the variance for each datapoint in :data:`y`.

   The covariance matrix for the parameters (:data:`c0`, :data:`c1`) is
   computed using the weights and returned via the parameters
   (:data:`cov00`, :data:`cov01`, :data:`cov11`).  The weighted sum of squares
   of the residuals from the best-fit line, :math:`\chi^2`, is returned in
   :data:`chisq`.

.. function:: int gsl_fit_linear_est (double x, double c0, double c1, double cov00, double cov01, double cov11, double * y, double * y_err)

   This function uses the best-fit linear regression coefficients
   :data:`c0`, :data:`c1` and their covariance
   :data:`cov00`, :data:`cov01`, :data:`cov11` to compute the fitted function
   :data:`y` and its standard deviation :data:`y_err` for the model :math:`Y = c_0 + c_1 X`
   at the point :data:`x`.

Linear regression without a constant term
-----------------------------------------

The functions described in this section can be used to perform
least-squares fits to a straight line model without a constant term,
:math:`Y = c_1 X`.

.. function:: int gsl_fit_mul (const double * x, const size_t xstride, const double * y, const size_t ystride, size_t n, double * c1, double * cov11, double * sumsq)

   This function computes the best-fit linear regression coefficient
   :data:`c1` of the model :math:`Y = c_1 X` for the datasets (:data:`x`,
   :data:`y`), two vectors of length :data:`n` with strides :data:`xstride` and
   :data:`ystride`.  The errors on :data:`y` are assumed unknown so the 
   variance of the parameter :data:`c1` is estimated from
   the scatter of the points around the best-fit line and returned via the
   parameter :data:`cov11`.  The sum of squares of the residuals from the
   best-fit line is returned in :data:`sumsq`.

.. function:: int gsl_fit_wmul (const double * x, const size_t xstride, const double * w, const size_t wstride, const double * y, const size_t ystride, size_t n, double * c1, double * cov11, double * sumsq)

   This function computes the best-fit linear regression coefficient
   :data:`c1` of the model :math:`Y = c_1 X` for the weighted datasets
   (:data:`x`, :data:`y`), two vectors of length :data:`n` with strides
   :data:`xstride` and :data:`ystride`.  The vector :data:`w`, of length :data:`n`
   and stride :data:`wstride`, specifies the weight of each datapoint. The
   weight is the reciprocal of the variance for each datapoint in :data:`y`.

   The variance of the parameter :data:`c1` is computed using the weights
   and returned via the parameter :data:`cov11`.  The weighted sum of
   squares of the residuals from the best-fit line, :math:`\chi^2`, is
   returned in :data:`chisq`.

.. function:: int gsl_fit_mul_est (double x, double c1, double cov11, double * y, double * y_err)

   This function uses the best-fit linear regression coefficient :data:`c1`
   and its covariance :data:`cov11` to compute the fitted function
   :data:`y` and its standard deviation :data:`y_err` for the model :math:`Y = c_1 X`
   at the point :data:`x`.

.. index::
   single: multi-parameter regression
   single: fits, multi-parameter linear

Multi-parameter regression
==========================

This section describes routines which perform least squares fits
to a linear model by minimizing the cost function

.. math:: \chi^2 = \sum_i w_i (y_i - \sum_j X_{ij} c_j)^2 = || y - Xc ||_W^2

where :math:`y` is a vector of :math:`n` observations, :math:`X` is an
:math:`n`-by-:math:`p` matrix of predictor variables, :math:`c`
is a vector of the :math:`p` unknown best-fit parameters to be estimated,
and :math:`||r||_W^2 = r^T W r`.
The matrix :math:`W = \diag(w_1,w_2,...,w_n)`
defines the weights or uncertainties of the observation vector.

This formulation can be used for fits to any number of functions and/or
variables by preparing the :math:`n`-by-:math:`p` matrix :math:`X`
appropriately.  For example, to fit to a :math:`p`-th order polynomial in
:data:`x`, use the following matrix,

.. math:: X_{ij} = x_i^j

where the index :math:`i` runs over the observations and the index
:math:`j` runs from 0 to :math:`p-1`.

To fit to a set of :math:`p` sinusoidal functions with fixed frequencies
:math:`\omega_1`, :math:`\omega_2`, :math:`\ldots`, :math:`\omega_p`, use,

.. math:: X_{ij} = \sin(\omega_j x_i)

To fit to :math:`p` independent variables :math:`x_1`, :math:`x_2`, :math:`\ldots`,
:math:`x_p`, use,

.. math:: X_{ij} = x_j(i)

where :math:`x_j(i)` is the :math:`i`-th value of the predictor variable
:math:`x_j`.

The solution of the general linear least-squares system requires an
additional working space for intermediate results, such as the singular
value decomposition of the matrix :math:`X`.

These functions are declared in the header file :file:`gsl_multifit.h`.

.. type:: gsl_multifit_linear_workspace

   This workspace contains internal variables for fitting multi-parameter models.

.. function:: gsl_multifit_linear_workspace * gsl_multifit_linear_alloc (const size_t n, const size_t p)

   This function allocates a workspace for fitting a model to a maximum of :data:`n`
   observations using a maximum of :data:`p` parameters. The user may later supply
   a smaller least squares system if desired. The size of the workspace is
   :math:`O(np + p^2)`.

.. function:: void gsl_multifit_linear_free (gsl_multifit_linear_workspace * work)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_multifit_linear_svd (const gsl_matrix * X, gsl_multifit_linear_workspace * work)

   This function performs a singular value decomposition of the
   matrix :data:`X` and stores the SVD factors internally in :data:`work`.

.. function:: int gsl_multifit_linear_bsvd (const gsl_matrix * X, gsl_multifit_linear_workspace * work)

   This function performs a singular value decomposition of the
   matrix :data:`X` and stores the SVD factors internally in :data:`work`.
   The matrix :data:`X` is first balanced by applying column scaling
   factors to improve the accuracy of the singular values.

.. function:: int gsl_multifit_linear (const gsl_matrix * X, const gsl_vector * y, gsl_vector * c, gsl_matrix * cov, double * chisq, gsl_multifit_linear_workspace * work)

   This function computes the best-fit parameters :data:`c` of the model
   :math:`y = X c` for the observations :data:`y` and the matrix of
   predictor variables :data:`X`, using the preallocated workspace provided
   in :data:`work`.  The :math:`p`-by-:math:`p` variance-covariance matrix of the model parameters
   :data:`cov` is set to :math:`\sigma^2 (X^T X)^{-1}`, where :math:`\sigma` is
   the standard deviation of the fit residuals.
   The sum of squares of the residuals from the best-fit,
   :math:`\chi^2`, is returned in :data:`chisq`. If the coefficient of
   determination is desired, it can be computed from the expression
   :math:`R^2 = 1 - \chi^2 / TSS`, where the total sum of squares (TSS) of
   the observations :data:`y` may be computed from :func:`gsl_stats_tss`.

   The best-fit is found by singular value decomposition of the matrix
   :data:`X` using the modified Golub-Reinsch SVD algorithm, with column
   scaling to improve the accuracy of the singular values. Any components
   which have zero singular value (to machine precision) are discarded
   from the fit.

.. function:: int gsl_multifit_linear_tsvd (const gsl_matrix * X, const gsl_vector * y, const double tol, gsl_vector * c, gsl_matrix * cov, double * chisq, size_t * rank, gsl_multifit_linear_workspace * work)

   This function computes the best-fit parameters :data:`c` of the model
   :math:`y = X c` for the observations :data:`y` and the matrix of
   predictor variables :data:`X`, using a truncated SVD expansion.
   Singular values which satisfy :math:`s_i \le tol \times s_0`
   are discarded from the fit, where :math:`s_0` is the largest singular value.
   The :math:`p`-by-:math:`p` variance-covariance matrix of the model parameters
   :data:`cov` is set to :math:`\sigma^2 (X^T X)^{-1}`, where :math:`\sigma` is
   the standard deviation of the fit residuals.
   The sum of squares of the residuals from the best-fit,
   :math:`\chi^2`, is returned in :data:`chisq`. The effective rank
   (number of singular values used in solution) is returned in :data:`rank`.
   If the coefficient of
   determination is desired, it can be computed from the expression
   :math:`R^2 = 1 - \chi^2 / TSS`, where the total sum of squares (TSS) of
   the observations :data:`y` may be computed from :func:`gsl_stats_tss`.

.. function:: int gsl_multifit_wlinear (const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, gsl_vector * c, gsl_matrix * cov, double * chisq, gsl_multifit_linear_workspace * work)

   This function computes the best-fit parameters :data:`c` of the weighted
   model :math:`y = X c` for the observations :data:`y` with weights :data:`w`
   and the matrix of predictor variables :data:`X`, using the preallocated
   workspace provided in :data:`work`.  The :math:`p`-by-:math:`p` covariance matrix of the model
   parameters :data:`cov` is computed as :math:`(X^T W X)^{-1}`. The weighted
   sum of squares of the residuals from the best-fit, :math:`\chi^2`, is
   returned in :data:`chisq`. If the coefficient of determination is
   desired, it can be computed from the expression :math:`R^2 = 1 - \chi^2 / WTSS`,
   where the weighted total sum of squares (WTSS) of the
   observations :data:`y` may be computed from :func:`gsl_stats_wtss`.

.. function:: int gsl_multifit_wlinear_tsvd (const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, const double tol, gsl_vector * c, gsl_matrix * cov, double * chisq, size_t * rank, gsl_multifit_linear_workspace * work)

   This function computes the best-fit parameters :data:`c` of the weighted
   model :math:`y = X c` for the observations :data:`y` with weights :data:`w`
   and the matrix of predictor variables :data:`X`, using a truncated SVD expansion.
   Singular values which satisfy :math:`s_i \le tol \times s_0`
   are discarded from the fit, where :math:`s_0` is the largest singular value.
   The :math:`p`-by-:math:`p` covariance matrix of the model
   parameters :data:`cov` is computed as :math:`(X^T W X)^{-1}`. The weighted
   sum of squares of the residuals from the best-fit, :math:`\chi^2`, is
   returned in :data:`chisq`. The effective rank of the system (number of
   singular values used in the solution) is returned in :data:`rank`.
   If the coefficient of determination is
   desired, it can be computed from the expression :math:`R^2 = 1 - \chi^2 / WTSS`,
   where the weighted total sum of squares (WTSS) of the
   observations :data:`y` may be computed from :func:`gsl_stats_wtss`.

.. function:: int gsl_multifit_linear_est (const gsl_vector * x, const gsl_vector * c, const gsl_matrix * cov, double * y, double * y_err)

   This function uses the best-fit multilinear regression coefficients
   :data:`c` and their covariance matrix
   :data:`cov` to compute the fitted function value
   :data:`y` and its standard deviation :data:`y_err` for the model :math:`y = x.c` 
   at the point :data:`x`.

.. function:: int gsl_multifit_linear_residuals (const gsl_matrix * X, const gsl_vector * y, const gsl_vector * c, gsl_vector * r)

   This function computes the vector of residuals :math:`r = y - X c` for
   the observations :data:`y`, coefficients :data:`c` and matrix of predictor
   variables :data:`X`.

.. function:: size_t gsl_multifit_linear_rank (const double tol, const gsl_multifit_linear_workspace * work)

   This function returns the rank of the matrix :math:`X` which must first have its
   singular value decomposition computed. The rank is computed by counting the number
   of singular values :math:`\sigma_j` which satisfy :math:`\sigma_j > tol \times \sigma_0`,
   where :math:`\sigma_0` is the largest singular value.

.. index::
   single: ridge regression
   single: Tikhonov regression
   single: regression, ridge
   single: regression, Tikhonov
   single: least squares, regularized

.. _sec_regularized-regression:

Regularized regression
======================

Ordinary weighted least squares models seek a solution vector :math:`c`
which minimizes the residual

.. math:: \chi^2 = || y - Xc ||_W^2

where :math:`y` is the :math:`n`-by-:math:`1` observation vector,
:math:`X` is the :math:`n`-by-:math:`p` design matrix, :math:`c` is
the :math:`p`-by-:math:`1` solution vector,
:math:`W = \diag(w_1,...,w_n)` is the data weighting matrix,
and :math:`||r||_W^2 = r^T W r`.
In cases where the least squares matrix :math:`X` is ill-conditioned,
small perturbations (ie: noise) in the observation vector could lead to
widely different solution vectors :math:`c`.
One way of dealing with ill-conditioned matrices is to use a "truncated SVD"
in which small singular values, below some given tolerance, are discarded
from the solution. The truncated SVD method is available using the functions
:func:`gsl_multifit_linear_tsvd` and :func:`gsl_multifit_wlinear_tsvd`. Another way
to help solve ill-posed problems is to include a regularization term in the least squares
minimization

.. math:: \chi^2 = || y - Xc ||_W^2 + \lambda^2 || L c ||^2

for a suitably chosen regularization parameter :math:`\lambda` and
matrix :math:`L`. This type of regularization is known as Tikhonov, or ridge,
regression. In some applications, :math:`L` is chosen as the identity matrix, giving
preference to solution vectors :math:`c` with smaller norms.
Including this regularization term leads to the explicit "normal equations" solution

.. only:: not texinfo

   .. math:: c = \left( X^T W X + \lambda^2 L^T L \right)^{-1} X^T W y

.. only:: texinfo

   ::

      c = ( X^T W X + \lambda^2 L^T L )^-1 X^T W y

which reduces to the ordinary least squares solution when :math:`L = 0`.
In practice, it is often advantageous to transform a regularized least
squares system into the form

.. only:: not texinfo

   .. math:: \chi^2 = || \tilde{y} - \tilde{X} \tilde{c} ||^2 + \lambda^2 || \tilde{c} ||^2

.. only:: texinfo

   ::

      \chi^2 = || y~ - X~ c~ ||^2 + \lambda^2 || c~ ||^2

This is known as the Tikhonov "standard form" and has the normal equations solution

.. only:: not texinfo

   .. math:: \tilde{c} = \left( \tilde{X}^T \tilde{X} + \lambda^2 I \right)^{-1} \tilde{X}^T \tilde{y}

.. only:: texinfo

   ::

      \tilde{c} = ( \tilde{X}^T \tilde{X} + \lambda^2 I )^{-1} \tilde{X}^T \tilde{y}

For an :math:`m`-by-:math:`p` matrix :math:`L` which is full rank and has :math:`m >= p` (ie: :math:`L` is
square or has more rows than columns), we can calculate the "thin" QR decomposition of :math:`L`, and
note that :math:`||L c|| = ||R c||` since the :math:`Q` factor will not change the norm. Since
:math:`R` is :math:`p`-by-:math:`p`, we can then use the transformation

.. only:: not texinfo

   .. math::

      \tilde{X} &= W^{1 \over 2} X R^{-1} \\
      \tilde{y} &= W^{1 \over 2} y \\
      \tilde{c} &= R c

.. only:: texinfo

   ::

      X~ = sqrt(W) X R^-1
      y~ = sqrt(W) y
      c~ = R c

to achieve the standard form. For a rectangular matrix :math:`L` with :math:`m < p`,
a more sophisticated approach is needed (see Hansen 1998, chapter 2.3).
In practice, the normal equations solution above is not desirable due to
numerical instabilities, and so the system is solved using the
singular value decomposition of the matrix :math:`\tilde{X}`.
The matrix :math:`L` is often chosen as the identity matrix, or as a first
or second finite difference operator, to ensure a smoothly varying
coefficient vector :math:`c`, or as a diagonal matrix to selectively damp
each model parameter differently. If :math:`L \ne I`, the user must first
convert the least squares problem to standard form using
:func:`gsl_multifit_linear_stdform1` or :func:`gsl_multifit_linear_stdform2`,
solve the system, and then backtransform the solution vector to recover
the solution of the original problem (see
:func:`gsl_multifit_linear_genform1` and :func:`gsl_multifit_linear_genform2`).

In many regularization problems, care must be taken when choosing
the regularization parameter :math:`\lambda`. Since both the
residual norm :math:`||y - X c||` and solution norm :math:`||L c||`
are being minimized, the parameter :math:`\lambda` represents
a tradeoff between minimizing either the residuals or the
solution vector. A common tool for visualizing the comprimise between
the minimization of these two quantities is known as the L-curve.
The L-curve is a log-log plot of the residual norm :math:`||y - X c||`
on the horizontal axis and the solution norm :math:`||L c||` on the
vertical axis. This curve nearly always as an :math:`L` shaped
appearance, with a distinct corner separating the horizontal
and vertical sections of the curve. The regularization parameter
corresponding to this corner is often chosen as the optimal
value. GSL provides routines to calculate the L-curve for all
relevant regularization parameters as well as locating the corner.

Another method of choosing the regularization parameter is known
as Generalized Cross Validation (GCV). This method is based on
the idea that if an arbitrary element :math:`y_i` is left out of the
right hand side, the resulting regularized solution should predict this element
accurately. This leads to choosing the parameter :math:`\lambda`
which minimizes the GCV function

.. only:: not texinfo

   .. math:: G(\lambda) = {||y - X c_{\lambda}||^2 \over \textrm{Tr}(I_n - X X_{\lambda}^I)^2}

.. only:: texinfo

   ::

      G(\lambda) = (||y - X c_{\lambda}||^2) / Tr(I_n - X X^I)^2

where :math:`X_{\lambda}^I` is the matrix which relates the solution :math:`c_{\lambda}`
to the right hand side :math:`y`, ie: :math:`c_{\lambda} = X_{\lambda}^I y`. GSL
provides routines to compute the GCV curve and its minimum.

For most applications, the steps required to solve a regularized least
squares problem are as follows:

#. Construct the least squares system (:math:`X`, :math:`y`, :math:`W`, :math:`L`)

#. Transform the system to standard form (:math:`\tilde{X}`, :math:`\tilde{y}`). This
   step can be skipped if :math:`L = I_p` and :math:`W = I_n`.

#. Calculate the SVD of :math:`\tilde{X}`.

#. Determine an appropriate regularization parameter :math:`\lambda` (using for example
   L-curve or GCV analysis).

#. Solve the standard form system using the chosen :math:`\lambda` and the SVD of :math:`\tilde{X}`.

#. Backtransform the standard form solution :math:`\tilde{c}` to recover the
   original solution vector :math:`c`.

.. function:: int gsl_multifit_linear_stdform1 (const gsl_vector * L, const gsl_matrix * X, const gsl_vector * y, gsl_matrix * Xs, gsl_vector * ys, gsl_multifit_linear_workspace * work)
              int gsl_multifit_linear_wstdform1 (const gsl_vector * L, const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, gsl_matrix * Xs, gsl_vector * ys, gsl_multifit_linear_workspace * work)

   These functions define a regularization matrix
   :math:`L = \diag(l_0,l_1,...,l_{p-1})`.
   The diagonal matrix element :math:`l_i` is provided by the
   :math:`i`-th element of the input vector :data:`L`.
   The :math:`n`-by-:math:`p` least squares matrix :data:`X` and
   vector :data:`y` of length :math:`n` are then
   converted to standard form as described above and the parameters
   (:math:`\tilde{X}`, :math:`\tilde{y}`) are stored in :data:`Xs` and :data:`ys`
   on output.  :data:`Xs` and :data:`ys` have the same dimensions as
   :data:`X` and :data:`y`. Optional data weights may be supplied in the
   vector :data:`w` of length :math:`n`. In order to apply this transformation,
   :math:`L^{-1}` must exist and so none of the :math:`l_i`
   may be zero. After the standard form system has been solved,
   use :func:`gsl_multifit_linear_genform1` to recover the original solution vector.
   It is allowed to have :data:`X` = :data:`Xs` and :data:`y` = :data:`ys` for an in-place transform.
   In order to perform a weighted regularized fit with :math:`L = I`, the user may
   call :func:`gsl_multifit_linear_applyW` to convert to standard form.

.. function:: int gsl_multifit_linear_L_decomp (gsl_matrix * L, gsl_vector * tau)

   This function factors the :math:`m`-by-:math:`p` regularization matrix
   :data:`L` into a form needed for the later transformation to standard form. :data:`L`
   may have any number of rows :math:`m`. If :math:`m \ge p` the QR decomposition of
   :data:`L` is computed and stored in :data:`L` on output. If :math:`m < p`, the QR decomposition
   of :math:`L^T` is computed and stored in :data:`L` on output. On output,
   the Householder scalars are stored in the vector :data:`tau` of size :math:`MIN(m,p)`.
   These outputs will be used by :func:`gsl_multifit_linear_wstdform2` to complete the
   transformation to standard form.

.. function:: int gsl_multifit_linear_stdform2 (const gsl_matrix * LQR, const gsl_vector * Ltau, const gsl_matrix * X, const gsl_vector * y, gsl_matrix * Xs, gsl_vector * ys, gsl_matrix * M, gsl_multifit_linear_workspace * work)
              int gsl_multifit_linear_wstdform2 (const gsl_matrix * LQR, const gsl_vector * Ltau, const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, gsl_matrix * Xs, gsl_vector * ys, gsl_matrix * M, gsl_multifit_linear_workspace * work)

   These functions convert the least squares system (:data:`X`, :data:`y`, :data:`W`, :math:`L`) to standard
   form (:math:`\tilde{X}`, :math:`\tilde{y}`) which are stored in :data:`Xs` and :data:`ys`
   respectively. The :math:`m`-by-:math:`p` regularization matrix :data:`L` is specified by the inputs
   :data:`LQR` and :data:`Ltau`, which are outputs from :func:`gsl_multifit_linear_L_decomp`.
   The dimensions of the standard form parameters (:math:`\tilde{X}`, :math:`\tilde{y}`)
   depend on whether :math:`m` is larger or less than :math:`p`. For :math:`m \ge p`,
   :data:`Xs` is :math:`n`-by-:math:`p`, :data:`ys` is :math:`n`-by-1, and :data:`M` is
   not used. For :math:`m < p`, :data:`Xs` is :math:`(n - p + m)`-by-:math:`m`,
   :data:`ys` is :math:`(n - p + m)`-by-1, and :data:`M` is additional :math:`n`-by-:math:`p` workspace,
   which is required to recover the original solution vector after the system has been
   solved (see :func:`gsl_multifit_linear_genform2`). Optional data weights may be supplied in the
   vector :data:`w` of length :math:`n`, where :math:`W = \diag(w)`.

.. function:: int gsl_multifit_linear_solve (const double lambda, const gsl_matrix * Xs, const gsl_vector * ys, gsl_vector * cs, double * rnorm, double * snorm, gsl_multifit_linear_workspace * work)

   This function computes the regularized best-fit parameters :math:`\tilde{c}`
   which minimize the cost function
   :math:`\chi^2 = || \tilde{y} - \tilde{X} \tilde{c} ||^2 + \lambda^2 || \tilde{c} ||^2`
   which is in standard form. The least squares system must therefore be converted
   to standard form prior to calling this function.
   The observation vector :math:`\tilde{y}` is provided in :data:`ys` and the matrix of
   predictor variables :math:`\tilde{X}` in :data:`Xs`. The solution vector :math:`\tilde{c}` is
   returned in :data:`cs`, which has length min(:math:`m,p`). The SVD of :data:`Xs` must be computed prior
   to calling this function, using :func:`gsl_multifit_linear_svd`.
   The regularization parameter :math:`\lambda` is provided in :data:`lambda`.
   The residual norm :math:`|| \tilde{y} - \tilde{X} \tilde{c} || = ||y - X c||_W`
   is returned in :data:`rnorm`.
   The solution norm :math:`|| \tilde{c} || = ||L c||` is returned in
   :data:`snorm`.

.. function:: int gsl_multifit_linear_genform1 (const gsl_vector * L, const gsl_vector * cs, gsl_vector * c, gsl_multifit_linear_workspace * work)

   After a regularized system has been solved with
   :math:`L = \diag(\l_0,\l_1,...,\l_{p-1})`,
   this function backtransforms the standard form solution vector :data:`cs`
   to recover the solution vector of the original problem :data:`c`. The
   diagonal matrix elements :math:`l_i` are provided in
   the vector :data:`L`. It is allowed to have :data:`c` = :data:`cs` for an
   in-place transform.

.. function:: int gsl_multifit_linear_genform2 (const gsl_matrix * LQR, const gsl_vector * Ltau, const gsl_matrix * X, const gsl_vector * y, const gsl_vector * cs, const gsl_matrix * M, gsl_vector * c, gsl_multifit_linear_workspace * work)
              int gsl_multifit_linear_wgenform2 (const gsl_matrix * LQR, const gsl_vector * Ltau, const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, const gsl_vector * cs, const gsl_matrix * M, gsl_vector * c, gsl_multifit_linear_workspace * work)

   After a regularized system has been solved with a general rectangular matrix :math:`L`,
   specified by (:data:`LQR`, :data:`Ltau`), this function backtransforms the standard form solution :data:`cs`
   to recover the solution vector of the original problem, which is stored in :data:`c`,
   of length :math:`p`. The original least squares matrix and observation vector are provided in
   :data:`X` and :data:`y` respectively. :data:`M` is the matrix computed by
   :func:`gsl_multifit_linear_stdform2`. For weighted fits, the weight vector
   :data:`w` must also be supplied.

.. function:: int gsl_multifit_linear_applyW (const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, gsl_matrix * WX, gsl_vector * Wy)

   For weighted least squares systems with :math:`L = I`, this function may be used to
   convert the system to standard form by applying the weight matrix :math:`W = \diag(w)`
   to the least squares matrix :data:`X` and observation vector :data:`y`. On output, :data:`WX`
   is equal to :math:`W^{1/2} X` and :data:`Wy` is equal to :math:`W^{1/2} y`. It is allowed
   for :data:`WX` = :data:`X` and :data:`Wy` = :data:`y` for an in-place transform.

.. function:: int gsl_multifit_linear_lreg (const double smin, const double smax, gsl_vector * reg_param)

   This function computes a set of possible regularization parameters for L-curve analysis
   and stores them in the output vector :data:`reg_param` of length :math:`k`, where :math:`k` is
   the number of desired points on the L-curve. The regularization parameters are equally logarithmically
   distributed between the provided values :data:`smin` and :data:`smax`, which typically correspond to
   the minimum and maximum singular value of the least squares matrix in standard form.
   The regularization parameters are calculated as,

   .. math:: \lambda_i = S_{min} \left( \frac{S_{max}}{S_{min}} \right)^{\frac{i-1}{k-1}}, \quad i = 1, \dots, k

.. function:: int gsl_multifit_linear_lcurve (const gsl_vector * y, gsl_vector * reg_param, gsl_vector * rho, gsl_vector * eta, gsl_multifit_linear_workspace * work)

   This function computes the L-curve for a least squares system
   using the right hand side vector :data:`y` and the SVD decomposition
   of the least squares matrix :data:`X`, which must be provided
   to :func:`gsl_multifit_linear_svd` prior to
   calling this function. The output vectors :data:`reg_param`,
   :data:`rho`, and :data:`eta` must all be the same size, and will
   contain the regularization parameters :math:`\lambda_i`, residual norms
   :math:`||y - X c_i||`, and solution norms :math:`|| L c_i ||`
   which compose the L-curve, where :math:`c_i` is the regularized
   solution vector corresponding to :math:`\lambda_i`.
   The user may determine the number of points on the L-curve by
   adjusting the size of these input arrays. The regularization
   parameters :math:`\lambda_i` are estimated from the singular values
   of :data:`X`, and chosen to represent the most relevant portion of
   the L-curve.

.. function:: int gsl_multifit_linear_lcurvature (const gsl_vector * y, const gsl_vector * reg_param, const gsl_vector * rho, const gsl_vector * eta, gsl_vector * kappa, gsl_multifit_linear_workspace * work)

   This function computes the curvature of the L-curve
   :math:`(\log{\hat{\rho}(\lambda)}, \log{\hat{\eta}(\lambda)})`, where
   :math:`\hat{\rho}(\lambda) = \log{||y - X c_{\lambda}||}` and
   :math:`\hat{\eta}(\lambda) = \log{|| L c_{\lambda} ||}`.
   This function uses the right hand side vector :data:`y`,
   the vector of regularization parameters, :data:`reg_param`,
   vector of residual norms, :data:`rho`, and vector of solution norms, :data:`eta`.
   The arrays :data:`reg_param`, :data:`rho`, and :data:`eta` can be computed by
   :func:`gsl_multifit_linear_lcurve`.  The curvature is defined as

   .. math:: \kappa(\lambda) = \frac{\hat{\rho}' \hat{\eta}'' - \hat{\rho}'' \hat{\eta}'}{\left( (\hat{\rho}')^2 + (\hat{\eta}')^2 \right)^{\frac{3}{2}}}

   The curvature values are stored in :data:`kappa` on output. The function
   :func:`gsl_multifit_linear_svd` must be called on the least squares matrix
   :math:`X` prior to calling this function.

.. function:: int gsl_multifit_linear_lcurvature_menger(const gsl_vector * rho, const gsl_vector * eta, gsl_vector * kappa)

   This function computes the Menger curvature of the L-curve
   :math:`(\log{\hat{\rho}(\lambda)}, \log{\hat{\eta}(\lambda)})`, where
   :math:`\hat{\rho}(\lambda) = \log{||y - X c_{\lambda}||}` and
   :math:`\hat{\eta}(\lambda) = \log{|| L c_{\lambda} ||}`.
   The function computes the Menger curvature for each consecutive triplet of points,
   :math:`\kappa_i = 1/R_i`, where :math:`R_i` is the radius of the unique
   circle fitted to the points
   :math:`(\hat{\rho}_{i-1}, \hat{\eta}_{i-1}), (\hat{\rho}_i, \hat{\eta}_i), (\hat{\rho}_{i+1}, \hat{\eta}_{i+1})`.

   The vector of residual norms :math:`|| y - X c_{\lambda} ||` is provided in :data:`rho`,
   the vector of solution norms :math:`|| L c_{\lambda} ||` is provided in :data:`eta`.
   These arrays can be calculated by calling :func:`gsl_multifit_linear_lcurve`.
   The Menger curvature output is stored in :data:`kappa`. The Menger curvature is
   an approximation to the curvature calculated by :func:`gsl_multifit_linear_lcurvature` but
   may be faster to calculate.

.. function:: int gsl_multifit_linear_lcorner (const gsl_vector * rho, const gsl_vector * eta, size_t * idx)

   This function attempts to locate the corner of the L-curve
   :math:`(||y - X c||, ||L c||)` defined by the :data:`rho` and :data:`eta`
   input arrays respectively. The corner is defined as the point of maximum
   curvature of the L-curve in log-log scale. The :data:`rho` and :data:`eta`
   arrays can be outputs of :func:`gsl_multifit_linear_lcurve`. The
   algorithm used simply fits a circle to 3 consecutive points on the L-curve
   and uses the circle's radius to determine the curvature at
   the middle point. Therefore, the input array sizes must be
   :math:`\ge 3`. With more points provided for the L-curve, a better
   estimate of the curvature can be obtained. The array index
   corresponding to maximum curvature (ie: the corner) is returned
   in :data:`idx`. If the input arrays contain colinear points,
   this function could fail and return :macro:`GSL_EINVAL`.

.. function:: int gsl_multifit_linear_lcorner2 (const gsl_vector * reg_param, const gsl_vector * eta, size_t * idx)

   This function attempts to locate the corner of an alternate L-curve
   :math:`(\lambda^2, ||L c||^2)` studied by Rezghi and Hosseini, 2009.
   This alternate L-curve can provide better estimates of the
   regularization parameter for smooth solution vectors. The regularization
   parameters :math:`\lambda` and solution norms :math:`||L c||` are provided
   in the :data:`reg_param` and :data:`eta` input arrays respectively. The
   corner is defined as the point of maximum curvature of this
   alternate L-curve in linear scale. The :data:`reg_param` and :data:`eta`
   arrays can be outputs of :func:`gsl_multifit_linear_lcurve`. The
   algorithm used simply fits a circle to 3 consecutive points on the L-curve
   and uses the circle's radius to determine the curvature at
   the middle point. Therefore, the input array sizes must be
   :math:`\ge 3`. With more points provided for the L-curve, a better
   estimate of the curvature can be obtained. The array index
   corresponding to maximum curvature (ie: the corner) is returned
   in :data:`idx`. If the input arrays contain colinear points,
   this function could fail and return :macro:`GSL_EINVAL`.

.. function:: int gsl_multifit_linear_gcv_init(const gsl_vector * y, gsl_vector * reg_param, gsl_vector * UTy, double * delta0, gsl_multifit_linear_workspace * work)

   This function performs some initialization in preparation for computing
   the GCV curve and its minimum. The right hand side vector is provided
   in :data:`y`. On output, :data:`reg_param` is set to a vector of regularization
   parameters in decreasing order and may be of any size. The vector
   :data:`UTy` of size :math:`p` is set to :math:`U^T y`. The parameter
   :data:`delta0` is needed for subsequent steps of the GCV calculation.

.. function:: int gsl_multifit_linear_gcv_curve(const gsl_vector * reg_param, const gsl_vector * UTy, const double delta0, gsl_vector * G, gsl_multifit_linear_workspace * work)

   This funtion calculates the GCV curve :math:`G(\lambda)` and stores it in
   :data:`G` on output, which must be the same size as :data:`reg_param`. The
   inputs :data:`reg_param`, :data:`UTy` and :data:`delta0` are computed in
   :func:`gsl_multifit_linear_gcv_init`.

.. function:: int gsl_multifit_linear_gcv_min(const gsl_vector * reg_param, const gsl_vector * UTy, const gsl_vector * G, const double delta0, double * lambda, gsl_multifit_linear_workspace * work)

   This function computes the value of the regularization parameter
   which minimizes the GCV curve :math:`G(\lambda)` and stores it in
   :data:`lambda`. The input :data:`G` is calculated by
   :func:`gsl_multifit_linear_gcv_curve` and the inputs
   :data:`reg_param`, :data:`UTy` and :data:`delta0` are computed by
   :func:`gsl_multifit_linear_gcv_init`.

.. function:: double gsl_multifit_linear_gcv_calc(const double lambda, const gsl_vector * UTy, const double delta0, gsl_multifit_linear_workspace * work)

   This function returns the value of the GCV curve :math:`G(\lambda)` corresponding
   to the input :data:`lambda`.

.. function:: int gsl_multifit_linear_gcv(const gsl_vector * y, gsl_vector * reg_param, gsl_vector * G, double * lambda, double * G_lambda, gsl_multifit_linear_workspace * work)

   This function combines the steps :code:`gcv_init`, :code:`gcv_curve`,
   and :code:`gcv_min` defined above into a single function. The input
   :data:`y` is the right hand side vector. On output, :data:`reg_param` and
   :data:`G`, which must be the same size, are set to vectors of
   :math:`\lambda` and :math:`G(\lambda)` values respectively. The
   output :data:`lambda` is set to the optimal value of :math:`\lambda`
   which minimizes the GCV curve. The minimum value of the GCV curve is
   returned in :data:`G_lambda`.

.. function:: int gsl_multifit_linear_Lk (const size_t p, const size_t k, gsl_matrix * L)

   This function computes the discrete approximation to the derivative operator :math:`L_k` of
   order :data:`k` on a regular grid of :data:`p` points and stores it in :data:`L`. The dimensions of :data:`L` are
   :math:`(p-k)`-by-:math:`p`.

.. function:: int gsl_multifit_linear_Lsobolev (const size_t p, const size_t kmax, const gsl_vector * alpha, gsl_matrix * L, gsl_multifit_linear_workspace * work)

   This function computes the regularization matrix :data:`L` corresponding to the
   weighted Sobolov norm
   :math:`||L c||^2 = \sum_k \alpha_k^2 ||L_k c||^2` where :math:`L_k` approximates
   the derivative operator of order :math:`k`. This regularization norm can be useful
   in applications where it is necessary to smooth several derivatives of the solution.
   :data:`p` is the number of model parameters, :data:`kmax` is the highest derivative
   to include in the summation above, and :data:`alpha` is the vector of weights of
   size :data:`kmax` + 1, where :code:`alpha[k]` = :math:`\alpha_k` is the weight
   assigned to the derivative of order :math:`k`.  The output matrix :data:`L` is size
   :data:`p`-by-:data:`p` and upper triangular.

.. function:: double gsl_multifit_linear_rcond (const gsl_multifit_linear_workspace * work)

   This function returns the reciprocal condition number of the least squares matrix :math:`X`,
   defined as the ratio of the smallest and largest singular values,
   rcond = :math:`\sigma_{min}/\sigma_{max}`.
   The routine :func:`gsl_multifit_linear_svd` must first be called to compute the SVD
   of :math:`X`.

.. index::
   single: robust regression
   single: regression, robust
   single: least squares, robust

Robust linear regression
========================

Ordinary least squares (OLS) models are often heavily influenced by the presence of outliers.
Outliers are data points which do not follow the general trend of the other observations,
although there is strictly no precise definition of an outlier. Robust linear regression
refers to regression algorithms which are robust to outliers. The most common type of
robust regression is M-estimation. The general M-estimator minimizes the objective function

.. math:: \sum_i \rho(e_i) = \sum_i \rho (y_i - Y(c, x_i))

where :math:`e_i = y_i - Y(c, x_i)` is the residual of the ith data point, and
:math:`\rho(e_i)` is a function which should have the following properties:

* :math:`\rho(e) \ge 0`
* :math:`\rho(0) = 0`
* :math:`\rho(-e) = \rho(e)`
* :math:`\rho(e_1) > \rho(e_2)` for :math:`|e_1| > |e_2|`

The special case of ordinary least squares is given by :math:`\rho(e_i) = e_i^2`.
Letting :math:`\psi = \rho'` be the derivative of :math:`\rho`, differentiating
the objective function with respect to the coefficients :math:`c`
and setting the partial derivatives to zero produces the system of equations

.. math:: \sum_i \psi(e_i) X_i = 0

where :math:`X_i` is a vector containing row :math:`i` of the design matrix :math:`X`.
Next, we define a weight function :math:`w(e) = \psi(e)/e`, and let
:math:`w_i = w(e_i)`:

.. math:: \sum_i w_i e_i X_i = 0

This system of equations is equivalent to solving a weighted ordinary least squares
problem, minimizing :math:`\chi^2 = \sum_i w_i e_i^2`. The weights however, depend
on the residuals :math:`e_i`, which depend on the coefficients :math:`c`, which depend
on the weights. Therefore, an iterative solution is used, called Iteratively Reweighted
Least Squares (IRLS).

#. Compute initial estimates of the coefficients :math:`c^{(0)}` using ordinary least squares

#. For iteration :math:`k`, form the residuals :math:`e_i^{(k)} = (y_i - X_i c^{(k-1)})/(t \sigma^{(k)} \sqrt{1 - h_i})`,
   where :math:`t` is a tuning constant depending on the choice of :math:`\psi`, and :math:`h_i` are the
   statistical leverages (diagonal elements of the matrix :math:`X (X^T X)^{-1} X^T`). Including :math:`t`
   and :math:`h_i` in the residual calculation has been shown to improve the convergence of the method.
   The residual standard deviation is approximated as :math:`\sigma^{(k)} = MAD / 0.6745`, where MAD is the
   Median-Absolute-Deviation of the :math:`n-p` largest residuals from the previous iteration.

#. Compute new weights :math:`w_i^{(k)} = \psi(e_i^{(k)})/e_i^{(k)}`.

#. Compute new coefficients :math:`c^{(k)}` by solving the weighted least squares problem with
   weights :math:`w_i^{(k)}`.

#. Steps 2 through 4 are iterated until the coefficients converge or until some maximum iteration
   limit is reached. Coefficients are tested for convergence using the critera:

  .. only:: not texinfo

     .. math:: |c_i^{(k)} - c_i^{(k-1)}| \le \epsilon \times \hbox{max}(|c_i^{(k)}|, |c_i^{(k-1)}|)

  .. only:: texinfo

     ::

        |c_i^(k) - c_i^(k-1)| <= \epsilon * max(|c_i^(k)|, |c_i^(k-1)|)

  for all :math:`0 \le i < p` where :math:`\epsilon` is a small tolerance factor.

The key to this method lies in selecting the function :math:`\psi(e_i)` to assign
smaller weights to large residuals, and larger weights to smaller residuals. As
the iteration proceeds, outliers are assigned smaller and smaller weights, eventually
having very little or no effect on the fitted model.

.. type:: gsl_multifit_robust_workspace

   This workspace is used for robust least squares fitting.

.. function:: gsl_multifit_robust_workspace * gsl_multifit_robust_alloc (const gsl_multifit_robust_type * T, const size_t n, const size_t p)

   This function allocates a workspace for fitting a model to :data:`n`
   observations using :data:`p` parameters. The size of the workspace
   is :math:`O(np + p^2)`. The type :data:`T` specifies the
   function :math:`\psi` and can be selected from the following choices.

   .. type:: gsl_multifit_robust_type

      .. var:: gsl_multifit_robust_type * gsl_multifit_robust_default

         This specifies the :data:`gsl_multifit_robust_bisquare` type (see below) and is a good
         general purpose choice for robust regression.

      .. var:: gsl_multifit_robust_type * gsl_multifit_robust_bisquare

         This is Tukey's biweight (bisquare) function and is a good general purpose choice for
         robust regression. The weight function is given by

         .. only:: not texinfo

            .. math::

               w(e) =
               \left\{
                 \begin{array}{cc}
                   (1 - e^2)^2, & |e| \le 1 \\
                    0, & |e| > 1
                 \end{array}
               \right.

         .. only:: texinfo

            ::

               w(e) = { (1 - e^2)^2, |e| <= 1
                      {     0,       |e| > 1

         and the default tuning constant is :math:`t = 4.685`.

      .. var:: gsl_multifit_robust_type * gsl_multifit_robust_cauchy

         This is Cauchy's function, also known as the Lorentzian function.
         This function does not guarantee a unique solution,
         meaning different choices of the coefficient vector :data:`c`
         could minimize the objective function. Therefore this option should
         be used with care. The weight function is given by

         .. only:: not texinfo

            .. math:: w(e) = {1 \over 1 + e^2}

         .. only:: texinfo

            ::

               w(e) = 1 / (1 + e^2)

         and the default tuning constant is :math:`t = 2.385`.

      .. var:: gsl_multifit_robust_type * gsl_multifit_robust_fair

         This is the fair :math:`\rho` function, which guarantees a unique solution and
         has continuous derivatives to three orders. The weight function is given by

         .. only:: not texinfo

            .. math:: w(e) = {1 \over 1 + |e|}

         .. only:: texinfo

            ::

               w(e) = 1 / (1 + |e|)

         and the default tuning constant is :math:`t = 1.400`.

      .. var:: gsl_multifit_robust_type * gsl_multifit_robust_huber

         This specifies Huber's :math:`\rho` function, which is a parabola in the vicinity of zero and
         increases linearly for a given threshold :math:`|e| > t`. This function is also considered
         an excellent general purpose robust estimator, however, occasional difficulties can
         be encountered due to the discontinuous first derivative of the :math:`\psi` function.
         The weight function is given by

         .. only:: not texinfo

            .. math::

               w(e) =
               \left\{
                 \begin{array}{cc}
                   1, & |e| \le 1 \\
                   {1 \over |e|}, & |e| > 1
                 \end{array}
               \right.

         .. only:: texinfo

            ::

               w(e) = 1/max(1,|e|)

         and the default tuning constant is :math:`t = 1.345`.

      .. var:: gsl_multifit_robust_type * gsl_multifit_robust_ols

         This specifies the ordinary least squares solution, which can be useful for quickly
         checking the difference between the various robust and OLS solutions. The
         weight function is given by

         .. math:: w(e) = 1

         and the default tuning constant is :math:`t = 1`.

      .. var:: gsl_multifit_robust_type * gsl_multifit_robust_welsch

         This specifies the Welsch function which can perform well in cases where the
         residuals have an exponential distribution. The weight function is given by

         .. math:: w(e) = \exp{(-e^2)}

         and the default tuning constant is :math:`t = 2.985`.

.. function:: void gsl_multifit_robust_free (gsl_multifit_robust_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: const char * gsl_multifit_robust_name (const gsl_multifit_robust_workspace * w)

   This function returns the name of the robust type :data:`T` specified to
   :func:`gsl_multifit_robust_alloc`.

.. function:: int gsl_multifit_robust_tune (const double tune, gsl_multifit_robust_workspace * w)

   This function sets the tuning constant :math:`t` used to adjust the residuals at each
   iteration to :data:`tune`.  Decreasing the tuning constant increases the downweight
   assigned to large residuals, while increasing the tuning constant decreases the
   downweight assigned to large residuals.

.. function:: int gsl_multifit_robust_maxiter (const size_t maxiter, gsl_multifit_robust_workspace * w)

   This function sets the maximum number of iterations in the iteratively
   reweighted least squares algorithm to :data:`maxiter`. By default,
   this value is set to 100 by :func:`gsl_multifit_robust_alloc`.

.. function:: int gsl_multifit_robust_weights (const gsl_vector * r, gsl_vector * wts, gsl_multifit_robust_workspace * w)

   This function assigns weights to the vector :data:`wts` using the residual vector
   :data:`r` and previously specified weighting function. The output weights are given
   by :math:`wts_i = w(r_i / (t \sigma))`, where the weighting functions :math:`w` are
   detailed in :func:`gsl_multifit_robust_alloc`. :math:`\sigma` is an estimate of the
   residual standard deviation based on the Median-Absolute-Deviation and :math:`t`
   is the tuning constant. This function is useful if the user wishes to implement
   their own robust regression rather than using
   the supplied :func:`gsl_multifit_robust` routine below.

.. function:: int gsl_multifit_robust (const gsl_matrix * X, const gsl_vector * y, gsl_vector * c, gsl_matrix * cov, gsl_multifit_robust_workspace * w)

   This function computes the best-fit parameters :data:`c` of the model
   :math:`y = X c` for the observations :data:`y` and the matrix of
   predictor variables :data:`X`, attemping to reduce the influence
   of outliers using the algorithm outlined above.
   The :math:`p`-by-:math:`p` variance-covariance matrix of the model parameters
   :data:`cov` is estimated as :math:`\sigma^2 (X^T X)^{-1}`, where :math:`\sigma` is
   an approximation of the residual standard deviation using the theory of robust
   regression. Special care must be taken when estimating :math:`\sigma` and
   other statistics such as :math:`R^2`, and so these
   are computed internally and are available by calling the function
   :func:`gsl_multifit_robust_statistics`.

   If the coefficients do not converge within the maximum iteration
   limit, the function returns :macro:`GSL_EMAXITER`. In this case,
   the current estimates of the coefficients and covariance matrix
   are returned in :data:`c` and :data:`cov` and the internal fit statistics
   are computed with these estimates.

.. function:: int gsl_multifit_robust_est (const gsl_vector * x, const gsl_vector * c, const gsl_matrix * cov, double * y, double * y_err)

   This function uses the best-fit robust regression coefficients
   :data:`c` and their covariance matrix
   :data:`cov` to compute the fitted function value
   :data:`y` and its standard deviation :data:`y_err` for the model :math:`y = x \cdot c` 
   at the point :data:`x`.

.. function:: int gsl_multifit_robust_residuals (const gsl_matrix * X, const gsl_vector * y, const gsl_vector * c, gsl_vector * r, gsl_multifit_robust_workspace * w)

   This function computes the vector of studentized residuals
   :math:`r_i = {y_i - (X c)_i \over \sigma \sqrt{1 - h_i}}` for
   the observations :data:`y`, coefficients :data:`c` and matrix of predictor
   variables :data:`X`. The routine :func:`gsl_multifit_robust` must
   first be called to compute the statisical leverages :math:`h_i` of
   the matrix :data:`X` and residual standard deviation estimate :math:`\sigma`.

.. function:: gsl_multifit_robust_stats gsl_multifit_robust_statistics (const gsl_multifit_robust_workspace * w)

   This function returns a structure containing relevant statistics from a robust regression.
   The function :func:`gsl_multifit_robust` must be called first to perform the regression and
   calculate these statistics.  The returned :type:`gsl_multifit_robust_stats` structure
   contains the following fields.

   .. type:: gsl_multifit_robust_stats

      :code:`double sigma_ols`
      
         This contains the standard deviation of the residuals as computed from ordinary
         least squares (OLS).

      :code:`double sigma_mad`
      
         This contains an estimate of the standard deviation of the final residuals using
         the Median-Absolute-Deviation statistic

      :code:`double sigma_rob`
      
         This contains an estimate of the standard deviation of the final residuals
         from the theory of robust regression (see Street et al, 1988).

      :code:`double sigma`
      
         This contains an estimate of the standard deviation of the final residuals
         by attemping to reconcile :code:`sigma_rob` and :code:`sigma_ols`
         in a reasonable way.

      :code:`double Rsq`
      
         This contains the :math:`R^2` coefficient of determination statistic using
         the estimate :code:`sigma`.

      :code:`double adj_Rsq`
      
         This contains the adjusted :math:`R^2` coefficient of determination statistic
         using the estimate :code:`sigma`.

      :code:`double rmse`
      
         This contains the root mean squared error of the final residuals

      :code:`double sse`
      
         This contains the residual sum of squares taking into account the robust
         covariance matrix.

      :code:`size_t dof`
      
         This contains the number of degrees of freedom :math:`n - p`

      :code:`size_t numit`
      
         Upon successful convergence, this contains the number of iterations performed

      :code:`gsl_vector * weights`
      
         This contains the final weight vector of length :data:`n`

      :code:`gsl_vector * r`
      
         This contains the final residual vector of length :data:`n`, :math:`r = y - X c`

.. index::
   single: large dense linear least squares
   single: linear least squares, large

Large dense linear systems
==========================

This module is concerned with solving large dense least squares systems
:math:`X c = y` where the :math:`n`-by-:math:`p` matrix
:math:`X` has :math:`n >> p` (ie: many more rows than columns).
This type of matrix is called a "tall skinny" matrix, and for
some applications, it may not be possible to fit the
entire matrix in memory at once to use the standard SVD approach.
Therefore, the algorithms in this module are designed to allow
the user to construct smaller blocks of the matrix :math:`X` and
accumulate those blocks into the larger system one at a time. The
algorithms in this module never need to store the entire matrix
:math:`X` in memory. The large linear least squares routines
support data weights and Tikhonov regularization, and are
designed to minimize the residual

.. math:: \chi^2 = || y - Xc ||_W^2 + \lambda^2 || L c ||^2

where :math:`y` is the :math:`n`-by-:math:`1` observation vector,
:math:`X` is the :math:`n`-by-:math:`p` design matrix, :math:`c` is
the :math:`p`-by-:math:`1` solution vector,
:math:`W = \diag(w_1,...,w_n)` is the data weighting matrix,
:math:`L` is an :math:`m`-by-:math:`p` regularization matrix,
:math:`\lambda` is a regularization parameter,
and :math:`||r||_W^2 = r^T W r`. In the discussion which follows,
we will assume that the system has been converted into Tikhonov
standard form,

.. only:: not texinfo

   .. math:: \chi^2 = || \tilde{y} - \tilde{X} \tilde{c} ||^2 + \lambda^2 || \tilde{c} ||^2

.. only:: texinfo

   ::

      \chi^2 = || y~ - X~ c~ ||^2 + \lambda^2 || c~ ||^2

and we will drop the tilde characters from the various parameters.
For a discussion of the transformation to standard form,
see :ref:`sec_regularized-regression`.

The basic idea is to partition the matrix :math:`X` and observation
vector :math:`y` as

.. only:: not texinfo

   .. math::

      \left(
        \begin{array}{c}
          X_1 \\
          X_2 \\
          X_3 \\
          \vdots \\
          X_k
        \end{array}
      \right)
      c =
      \left(
        \begin{array}{c}
          y_1 \\
          y_2 \\
          y_3 \\
          \vdots \\
          y_k
        \end{array}
      \right)

.. only:: texinfo

   ::

      [ X_1 ] c = [ y_1 ]
      [ X_2 ]     [ y_2 ]
      [ X_3 ]     [ y_3 ]
      [ ... ]     [ ... ]
      [ X_k ]     [ y_k ]

into :math:`k` blocks, where each block (:math:`X_i,y_i`) may have
any number of rows, but each :math:`X_i` has :math:`p` columns.
The sections below describe the methods available for solving
this partitioned system. The functions are declared in
the header file :file:`gsl_multilarge.h`.

.. index::
   single: large linear least squares, normal equations

Normal Equations Approach
-------------------------

The normal equations approach to the large linear least squares
problem described above is popular due to its speed and simplicity.
Since the normal equations solution to the problem is given by

.. only:: not texinfo

   .. math:: c = \left( X^T X + \lambda^2 I \right)^{-1} X^T y

.. only:: texinfo

   ::

      c = ( X^T X + \lambda^2 I )^-1 X^T y

only the :math:`p`-by-:math:`p` matrix :math:`X^T X` and
:math:`p`-by-1 vector :math:`X^T y` need to be stored. Using
the partition scheme described above, these are given by

.. only:: not texinfo

   .. math::

      X^T X &= \sum_i X_i^T X_i \\
      X^T y &= \sum_i X_i^T y_i

.. only:: texinfo

   ::

      X^T X = \sum_i X_i^T X_i
      X^T y = \sum_i X_i^T y_i

Since the matrix :math:`X^T X` is symmetric, only half of it
needs to be calculated. Once all of the blocks :math:`(X_i,y_i)`
have been accumulated into the final :math:`X^T X` and :math:`X^T y`,
the system can be solved with a Cholesky factorization of the
:math:`X^T X` matrix. The :math:`X^T X` matrix is first transformed via
a diagonal scaling transformation to attempt to reduce its condition
number as much as possible to recover a more accurate solution vector.
The normal equations approach is the fastest method for solving the
large least squares problem, and is accurate for well-conditioned
matrices :math:`X`. However, for ill-conditioned matrices, as is often
the case for large systems, this method can suffer from numerical
instabilities (see Trefethen and Bau, 1997).  The number of operations
for this method is :math:`O(np^2 + {1 \over 3}p^3)`.

.. index::
   single: large linear least squares, TSQR

Tall Skinny QR (TSQR) Approach
------------------------------

An algorithm which has better numerical stability for ill-conditioned
problems is known as the Tall Skinny QR (TSQR) method. This method
is based on computing the QR decomposition of the least squares
matrix,

.. only:: not texinfo

   .. math:: X = Q \begin{pmatrix} R \\ 0 \end{pmatrix}
   
.. only:: texinfo

   ::

      X = Q [ R; 0 ]

where :math:`Q` is an :math:`n`-by-:math:`n` orthogonal matrix,
and :math:`R` is a :math:`p`-by-:math:`p` upper triangular matrix.
If we define,

.. only:: not texinfo

   .. math:: z = Q^T y = \begin{pmatrix} z_1 \\ z_2 \end{pmatrix}

.. only:: texinfo

   ::

      z = Q^T y = [ z_1 ; z_2 ]

where :math:`z_1` is :math:`p`-by-1 and :math:`z_2` is :math:`(n-p)`-by-1,
then the residual becomes

.. math::
   
   \chi^2 &= \left|\left| y - X c \right|\right|^2 + \lambda^2 ||c||^2 \\
          &= \left|\left| Q^T y - \begin{pmatrix} R \\ 0 \end{pmatrix} c \right|\right|^2 + \lambda^2 || c ||^2 \\
          &= \left|\left| \begin{pmatrix} z_1 \\ z_2 \end{pmatrix} - \begin{pmatrix} R \\ 0 \end{pmatrix} c \right|\right|^2 + \lambda^2 || c ||^2 \\
          &= \left|\left| z_1 - R c \right|\right|^2 + || z_2 ||^2 + \lambda^2 || c ||^2 \\
          &= \left|\left| \begin{pmatrix} z_1 \\ 0 \end{pmatrix} - \begin{pmatrix} R \\ \lambda I_p \end{pmatrix} c \right|\right|^2 + || z_2 ||^2

Since :math:`z_2` does not depend on :math:`c`, :math:`\chi^2` is minimized
by solving the least squares system,

.. only:: not texinfo

   .. math::

      \left(
        \begin{array}{c}
          R \\
          \lambda I_p
        \end{array}
      \right) c =
      \left(
        \begin{array}{c}
          z_1 \\
          0
        \end{array}
      \right)

.. only:: texinfo

   ::

      [ R ; \lambda I ] c = [ Q^T b ; 0 ]

The matrix on the left hand side is now a much
smaller :math:`2p`-by-:math:`p` matrix which can
be solved with a standard SVD approach. The
:math:`Q` matrix is large, however it does not need to be
explicitly constructed. The TSQR algorithm
computes only the :math:`p`-by-:math:`p` matrix
:math:`R`, the :math:`p`-by-1 vector :math:`z_1`,
and the norm :math:`||z_2||`,
and updates these quantities as new blocks
are added to the system. Each time a new block of rows
(:math:`X_i,y_i`) is added, the algorithm performs a QR decomposition
of the matrix

.. only:: not texinfo

   .. math::

      \left(
        \begin{array}{c}
          R_{i-1} \\
          X_i
        \end{array}
      \right)

.. only:: texinfo

   ::

      [ R_(i-1) ; X_i ]

where :math:`R_{i-1}` is the upper triangular
:math:`R` factor for the matrix

.. only:: not texinfo

   .. math::

      \left(
        \begin{array}{c}
          X_1 \\
          \vdots \\
          X_{i-1}
        \end{array}
      \right)

.. only:: texinfo

   ::

      [ X_1 ; ... ; X_(i-1) ]

This QR decomposition is done efficiently taking into account
the sparse structure of :math:`R_{i-1}`. See Demmel et al, 2008 for
more details on how this is accomplished. The number
of operations for this method is :math:`O(2np^2 - {2 \over 3}p^3)`.

.. index::
   single: large linear least squares, steps

Large Dense Linear Systems Solution Steps
-----------------------------------------

The typical steps required to solve large regularized linear least
squares problems are as follows:

#. Choose the regularization matrix :math:`L`.

#. Construct a block of rows of the least squares matrix, right
   hand side vector, and weight vector (:math:`X_i`, :math:`y_i`, :math:`w_i`).

#. Transform the block to standard form (:math:`\tilde{X_i}`, :math:`\tilde{y_i}`). This
   step can be skipped if :math:`L = I` and :math:`W = I`.

#. Accumulate the standard form block (:math:`\tilde{X_i}`, :math:`\tilde{y_i}`) into
   the system.

#. Repeat steps 2-4 until the entire matrix and right hand side vector have
   been accumulated.

#. Determine an appropriate regularization parameter :math:`\lambda` (using for example
   L-curve analysis).

#. Solve the standard form system using the chosen :math:`\lambda`.

#. Backtransform the standard form solution :math:`\tilde{c}` to recover the
   original solution vector :math:`c`.

.. index::
   single: large linear least squares, routines

Large Dense Linear Least Squares Routines
-----------------------------------------

.. type:: gsl_multilarge_linear_workspace

   This workspace contains parameters for solving large linear least squares problems.

.. function:: gsl_multilarge_linear_workspace * gsl_multilarge_linear_alloc (const gsl_multilarge_linear_type * T, const size_t p)

   This function allocates a workspace for solving large linear least squares
   systems. The least squares matrix :math:`X` has :data:`p` columns,
   but may have any number of rows.
   
   .. type:: gsl_multilarge_linear_type

      The parameter :data:`T` specifies
      the method to be used for solving the large least squares system
      and may be selected from the following choices

      .. var:: gsl_multilarge_linear_type * gsl_multilarge_linear_normal

         This specifies the normal equations approach for
         solving the least squares system. This method is suitable
         in cases where performance is critical and it is known that the
         least squares matrix :math:`X` is well conditioned. The size
         of this workspace is :math:`O(p^2)`.

      .. var:: gsl_multilarge_linear_type * gsl_multilarge_linear_tsqr

         This specifies the sequential Tall Skinny QR (TSQR) approach for
         solving the least squares system. This method is a good
         general purpose choice for large systems, but requires about
         twice as many operations as the normal equations method for
         :math:`n >> p`. The size of this workspace is :math:`O(p^2)`.

.. function:: void gsl_multilarge_linear_free (gsl_multilarge_linear_workspace * w)

   This function frees the memory associated with the
   workspace :data:`w`.

.. function:: const char * gsl_multilarge_linear_name (gsl_multilarge_linear_workspace * w)

   This function returns a string pointer to the name
   of the multilarge solver.

.. function:: int gsl_multilarge_linear_reset (gsl_multilarge_linear_workspace * w)

   This function resets the workspace :data:`w` so
   it can begin to accumulate a new least squares
   system.

.. function:: int gsl_multilarge_linear_stdform1 (const gsl_vector * L, const gsl_matrix * X, const gsl_vector * y, gsl_matrix * Xs, gsl_vector * ys, gsl_multilarge_linear_workspace * work)
              int gsl_multilarge_linear_wstdform1 (const gsl_vector * L, const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, gsl_matrix * Xs, gsl_vector * ys, gsl_multilarge_linear_workspace * work)

   These functions define a regularization matrix
   :math:`L = \diag(l_0,l_1,...,l_{p-1})`.
   The diagonal matrix element :math:`l_i` is provided by the
   :math:`i`-th element of the input vector :data:`L`.
   The block (:data:`X`, :data:`y`) is converted to standard form and
   the parameters (:math:`\tilde{X}`, :math:`\tilde{y}`) are stored in :data:`Xs`
   and :data:`ys` on output.  :data:`Xs` and :data:`ys` have the same dimensions as
   :data:`X` and :data:`y`. Optional data weights may be supplied in the
   vector :data:`w`. In order to apply this transformation,
   :math:`L^{-1}` must exist and so none of the :math:`l_i`
   may be zero. After the standard form system has been solved,
   use :func:`gsl_multilarge_linear_genform1` to recover the original solution vector.
   It is allowed to have :data:`X` = :data:`Xs` and :data:`y` = :data:`ys` for an
   in-place transform.

.. function:: int gsl_multilarge_linear_L_decomp (gsl_matrix * L, gsl_vector * tau)

   This function calculates the QR decomposition of the :math:`m`-by-:math:`p`
   regularization matrix :data:`L`. :data:`L` must have :math:`m \ge p`.  On output,
   the Householder scalars are stored in the vector :data:`tau` of size :math:`p`.
   These outputs will be used by :func:`gsl_multilarge_linear_wstdform2` to complete the
   transformation to standard form.

.. function:: int gsl_multilarge_linear_stdform2 (const gsl_matrix * LQR, const gsl_vector * Ltau, const gsl_matrix * X, const gsl_vector * y, gsl_matrix * Xs, gsl_vector * ys, gsl_multilarge_linear_workspace * work)
              int gsl_multilarge_linear_wstdform2 (const gsl_matrix * LQR, const gsl_vector * Ltau, const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, gsl_matrix * Xs, gsl_vector * ys, gsl_multilarge_linear_workspace * work)

   These functions convert a block of rows (:data:`X`, :data:`y`, :data:`w`) to standard
   form (:math:`\tilde{X}`, :math:`\tilde{y}`) which are stored in :data:`Xs` and :data:`ys`
   respectively. :data:`X`, :data:`y`, and :data:`w` must all have the same number of rows.
   The :math:`m`-by-:math:`p` regularization matrix :data:`L` is specified by the inputs
   :data:`LQR` and :data:`Ltau`, which are outputs from :func:`gsl_multilarge_linear_L_decomp`.
   :data:`Xs` and :data:`ys` have the same dimensions as :data:`X` and :data:`y`. After the
   standard form system has been solved, use :func:`gsl_multilarge_linear_genform2` to
   recover the original solution vector. Optional data weights may be supplied in the
   vector :data:`w`, where :math:`W = \diag(w)`.

.. function:: int gsl_multilarge_linear_accumulate (gsl_matrix * X, gsl_vector * y, gsl_multilarge_linear_workspace * w)

   This function accumulates the standard form block (:math:`X,y`) into the
   current least squares system. :data:`X` and :data:`y` have the same number
   of rows, which can be arbitrary.  :data:`X` must have :math:`p` columns.
   For the TSQR method, :data:`X` and :data:`y` are destroyed on output.
   For the normal equations method, they are both unchanged.

.. function:: int gsl_multilarge_linear_solve (const double lambda, gsl_vector * c, double * rnorm, double * snorm, gsl_multilarge_linear_workspace * w)

   After all blocks (:math:`X_i,y_i`) have been accumulated into
   the large least squares system, this function will compute
   the solution vector which is stored in :data:`c` on output.
   The regularization parameter :math:`\lambda` is provided in
   :data:`lambda`. On output, :data:`rnorm` contains the residual norm
   :math:`||y - X c||_W` and :data:`snorm` contains the solution
   norm :math:`||L c||`.

.. function:: int gsl_multilarge_linear_genform1 (const gsl_vector * L, const gsl_vector * cs, gsl_vector * c, gsl_multilarge_linear_workspace * work)

   After a regularized system has been solved with
   :math:`L = \diag(\l_0,\l_1,...,\l_{p-1})`,
   this function backtransforms the standard form solution vector :data:`cs`
   to recover the solution vector of the original problem :data:`c`. The
   diagonal matrix elements :math:`l_i` are provided in
   the vector :data:`L`. It is allowed to have :data:`c` = :data:`cs` for an
   in-place transform.

.. function:: int gsl_multilarge_linear_genform2 (const gsl_matrix * LQR, const gsl_vector * Ltau, const gsl_vector * cs, gsl_vector * c, gsl_multilarge_linear_workspace * work)

   After a regularized system has been solved with a regularization matrix :math:`L`,
   specified by (:data:`LQR`, :data:`Ltau`), this function backtransforms the standard form
   solution :data:`cs`
   to recover the solution vector of the original problem, which is stored in :data:`c`,
   of length :math:`p`.

.. function:: int gsl_multilarge_linear_lcurve (gsl_vector * reg_param, gsl_vector * rho, gsl_vector * eta, gsl_multilarge_linear_workspace * work)

   This function computes the L-curve for a large least squares system
   after it has been fully accumulated into the workspace :data:`work`.
   The output vectors :data:`reg_param`, :data:`rho`, and :data:`eta` must all
   be the same size, and will contain the regularization parameters
   :math:`\lambda_i`, residual norms :math:`||y - X c_i||`, and solution
   norms :math:`|| L c_i ||` which compose the L-curve, where :math:`c_i`
   is the regularized solution vector corresponding to :math:`\lambda_i`.
   The user may determine the number of points on the L-curve by
   adjusting the size of these input arrays. For the TSQR method,
   the regularization parameters :math:`\lambda_i` are estimated from the
   singular values of the triangular :math:`R` factor. For the normal
   equations method, they are estimated from the eigenvalues of the
   :math:`X^T X` matrix.

.. function:: const gsl_matrix * gsl_multilarge_linear_matrix_ptr (const gsl_multilarge_linear_workspace * work)

   For the normal equations method, this function returns a pointer to the :math:`X^T X` matrix
   of size :math:`p`-by-:math:`p`. For the TSQR method, this function returns a pointer to the
   upper triangular :math:`R` matrix of size :math:`p`-by-:math:`p`.

.. function:: const gsl_vector * gsl_multilarge_linear_rhs_ptr (const gsl_multilarge_linear_workspace * work)

   For the normal equations method, this function returns a pointer to the :math:`p`-by-1 right hand
   side vector :math:`X^T y`. For the TSQR method, this function returns a pointer to a vector of length
   :math:`p+1`. The first :math:`p` elements of this vector contain :math:`z_1`, while the last element
   contains :math:`||z_2||`.

.. function:: int gsl_multilarge_linear_rcond (double * rcond, gsl_multilarge_linear_workspace * work)

   This function computes the reciprocal condition number, stored in
   :data:`rcond`, of the least squares matrix after it has been accumulated
   into the workspace :data:`work`. For the TSQR algorithm, this is
   accomplished by calculating the SVD of the :math:`R` factor, which
   has the same singular values as the matrix :math:`X`. For the normal
   equations method, this is done by computing the eigenvalues of
   :math:`X^T X`, which could be inaccurate for ill-conditioned matrices
   :math:`X`.

.. index:: least squares troubleshooting

Troubleshooting
===============

When using models based on polynomials, care should be taken when constructing the design matrix
:math:`X`. If the :math:`x` values are large, then the matrix :math:`X` could be ill-conditioned
since its columns are powers of :math:`x`, leading to unstable least-squares solutions.
In this case it can often help to center and scale the :math:`x` values using the mean and standard deviation:

.. only:: not texinfo

   .. math:: x' = {x - \mu(x) \over \sigma(x)}

.. only:: texinfo

   ::

      x' = (x - mu)/sigma

and then construct the :math:`X` matrix using the transformed values :math:`x'`.

Examples
========

The example programs in this section demonstrate the various linear regression methods.

Simple Linear Regression Example
--------------------------------

The following program computes a least squares straight-line fit to a
simple dataset, and outputs the best-fit line and its
associated one standard-deviation error bars.

.. include:: examples/fitting.c
   :code:

The following commands extract the data from the output of the program
and display it using the GNU plotutils "graph" utility::

  $ ./demo > tmp
  $ more tmp
  # best fit: Y = -106.6 + 0.06 X
  # covariance matrix:
  # [ 39602, -19.9
  #   -19.9, 0.01]
  # chisq = 0.8

  $ for n in data fit hi lo ; 
     do 
       grep "^$n" tmp | cut -d: -f2 > $n ; 
     done
  $ graph -T X -X x -Y y -y 0 20 -m 0 -S 2 -Ie data 
       -S 0 -I a -m 1 fit -m 2 hi -m 2 lo

The result is shown in :numref:`fig_fit-wlinear`.

.. _fig_fit-wlinear:

.. figure:: /images/fit-wlinear.png

   Straight line fit with 1-:math:`\sigma` error bars

Multi-parameter Linear Regression Example
-----------------------------------------

The following program performs a quadratic fit :math:`y = c_0 + c_1 x + c_2 x^2`
to a weighted dataset using the generalised linear fitting function
:func:`gsl_multifit_wlinear`.  The model matrix :math:`X` for a quadratic
fit is given by,

.. only:: not texinfo

   .. math::

      X =
      \left(
        \begin{array}{ccc}
          1 & x_0 & x_0^2 \\
          1 & x_1 & x_1^2 \\
          1 & x_2 & x_2^2 \\
          \dots & \dots & \dots
        \end{array}
      \right)

.. only:: texinfo

   ::

      X = [ 1   , x_0  , x_0^2 ;
            1   , x_1  , x_1^2 ;
            1   , x_2  , x_2^2 ;
            ... , ...  , ...   ]

where the column of ones corresponds to the constant term :math:`c_0`.
The two remaining columns corresponds to the terms :math:`c_1 x` and
:math:`c_2 x^2`.

The program reads :data:`n` lines of data in the format (:data:`x`, :data:`y`,
:data:`err`) where :data:`err` is the error (standard deviation) in the
value :data:`y`.

.. include:: examples/fitting2.c
   :code:

A suitable set of data for fitting can be generated using the following
program.  It outputs a set of points with gaussian errors from the curve
:math:`y = e^x` in the region :math:`0 < x < 2`.

.. include:: examples/fitting3.c
   :code:

The data can be prepared by running the resulting executable program::

  $ GSL_RNG_TYPE=mt19937_1999 ./generate > exp.dat
  $ more exp.dat
  0.1 0.97935 0.110517
  0.2 1.3359 0.12214
  0.3 1.52573 0.134986
  0.4 1.60318 0.149182
  0.5 1.81731 0.164872
  0.6 1.92475 0.182212
  ....

To fit the data use the previous program, with the number of data points
given as the first argument.  In this case there are 19 data points::

  $ ./fit 19 < exp.dat
  0.1 0.97935 +/- 0.110517
  0.2 1.3359 +/- 0.12214
  ...
  # best fit: Y = 1.02318 + 0.956201 X + 0.876796 X^2
  # covariance matrix:
  [ +1.25612e-02, -3.64387e-02, +1.94389e-02  
    -3.64387e-02, +1.42339e-01, -8.48761e-02  
    +1.94389e-02, -8.48761e-02, +5.60243e-02 ]
  # chisq = 23.0987

The parameters of the quadratic fit match the coefficients of the
expansion of :math:`e^x`, taking into account the errors on the
parameters and the :math:`O(x^3)` difference between the exponential and
quadratic functions for the larger values of :math:`x`.  The errors on
the parameters are given by the square-root of the corresponding
diagonal elements of the covariance matrix.  The chi-squared per degree
of freedom is 1.4, indicating a reasonable fit to the data.

:numref:`fig_fit-wlinear2` shows the resulting fit.

.. _fig_fit-wlinear2:

.. figure:: /images/fit-wlinear2.png

   Weighted fit example with error bars

Regularized Linear Regression Example 1
---------------------------------------

The next program demonstrates the difference between ordinary and
regularized least squares when the design matrix is near-singular.
In this program, we generate two random normally distributed variables
:math:`u` and :math:`v`, with :math:`v = u + noise` so that :math:`u`
and :math:`v` are nearly colinear. We then set a third dependent
variable :math:`y = u + v + noise` and solve for the coefficients
:math:`c_1,c_2` of the model :math:`Y(c_1,c_2) = c_1 u + c_2 v`.
Since :math:`u \approx v`, the design matrix :math:`X` is nearly
singular, leading to unstable ordinary least squares solutions.

Here is the program output::

  matrix condition number = 1.025113e+04

  === Unregularized fit ===
  best fit: y = -43.6588 u + 45.6636 v
  residual norm = 31.6248
  solution norm = 63.1764
  chisq/dof = 1.00213

  === Regularized fit (L-curve) ===
  optimal lambda: 4.51103
  best fit: y = 1.00113 u + 1.0032 v
  residual norm = 31.6547
  solution norm = 1.41728
  chisq/dof = 1.04499

  === Regularized fit (GCV) ===
  optimal lambda: 0.0232029
  best fit: y = -19.8367 u + 21.8417 v
  residual norm = 31.6332
  solution norm = 29.5051
  chisq/dof = 1.00314

We see that the ordinary least squares solution is completely wrong,
while the L-curve regularized method with the optimal
:math:`\lambda = 4.51103` finds the correct solution
:math:`c_1 \approx c_2 \approx 1`. The GCV regularized method finds
a regularization parameter :math:`\lambda = 0.0232029` which is too
small to give an accurate solution, although it performs better than OLS.
The L-curve and its computed corner, as well as the GCV curve and its
minimum are plotted in :numref:`fig_regularized`.

.. _fig_regularized:

.. figure:: /images/regularized.png

   L-curve and GCV curve for example program.

The program is given below.

.. include:: examples/fitreg.c
   :code:

Regularized Linear Regression Example 2
---------------------------------------

The following example program minimizes the cost function

.. math:: ||y - X c||^2 + \lambda^2 ||x||^2

where :math:`X` is the :math:`10`-by-:math:`8` Hilbert matrix whose
entries are given by

.. only:: not texinfo

   .. math:: X_{ij} = {1 \over i + j - 1}

.. only:: texinfo

   ::

      X_{ij} = 1 / (i + j - 1)

and the right hand side vector is given by
:math:`y = [1,-1,1,-1,1,-1,1,-1,1,-1]^T`. Solutions
are computed for :math:`\lambda = 0` (unregularized) as
well as for optimal parameters :math:`\lambda` chosen by
analyzing the L-curve and GCV curve.

Here is the program output::

  matrix condition number = 3.565872e+09

  === Unregularized fit ===
  residual norm = 2.15376
  solution norm = 2.92217e+09
  chisq/dof = 2.31934

  === Regularized fit (L-curve) ===
  optimal lambda: 7.11407e-07
  residual norm = 2.60386
  solution norm = 424507
  chisq/dof = 3.43565

  === Regularized fit (GCV) ===
  optimal lambda: 1.72278
  residual norm = 3.1375
  solution norm = 0.139357
  chisq/dof = 4.95076

Here we see the unregularized solution results in a large solution
norm due to the ill-conditioned matrix. The L-curve solution finds
a small value of :math:`\lambda = 7.11e-7` which still results in
a badly conditioned system and a large solution norm. The GCV method
finds a parameter :math:`\lambda = 1.72` which results in a well-conditioned
system and small solution norm.

The L-curve and its computed corner, as well as the GCV curve and its
minimum are plotted in :numref:`fig_regularized2`.

.. _fig_regularized2:

.. figure:: /images/regularized2.png

   L-curve and GCV curve for example program.

The program is given below.

.. include:: examples/fitreg2.c
   :code:

Robust Linear Regression Example
--------------------------------

The next program demonstrates the advantage of robust least squares on
a dataset with outliers. The program generates linear :math:`(x,y)`
data pairs on the line :math:`y = 1.45 x + 3.88`, adds some random
noise, and inserts 3 outliers into the dataset. Both the robust
and ordinary least squares (OLS) coefficients are computed for
comparison.

.. include:: examples/robfit.c
   :code:

The output from the program is shown in :numref:`fig_robust`.

.. _fig_robust:

.. figure:: /images/robust.png

   Linear fit to dataset with outliers.

Large Dense Linear Regression Example
-------------------------------------

The following program demostrates the large dense linear least squares
solvers. This example is adapted from Trefethen and Bau,
and fits the function :math:`f(t) = \exp{(\sin^3{(10t)}})` on
the interval :math:`[0,1]` with a degree 15 polynomial. The
program generates :math:`n = 50000` equally spaced points
:math:`t_i` on this interval, calculates the function value
and adds random noise to determine the observation value
:math:`y_i`. The entries of the least squares matrix are
:math:`X_{ij} = t_i^j`, representing a polynomial fit. The
matrix is highly ill-conditioned, with a condition number
of about :math:`2.4 \cdot 10^{11}`. The program accumulates the
matrix into the least squares system in 5 blocks, each with
10000 rows. This way the full matrix :math:`X` is never
stored in memory. We solve the system with both the
normal equations and TSQR methods. The results are shown
in :numref:`fig_multilarge`. In the top left plot, the TSQR
solution provides a reasonable agreement to the exact solution,
while the normal equations method fails completely since the
Cholesky factorization fails due to the ill-conditioning of the matrix.
In the bottom left plot, we show the L-curve calculated from TSQR, which exhibits multiple corners.
In the top right panel, we plot a regularized solution using
:math:`\lambda = 10^{-5}`. The TSQR and normal solutions now agree,
however they are unable to provide a good fit due to the damping.
This indicates that for some ill-conditioned
problems, regularizing the normal equations does not improve the
solution. This is further illustrated in the bottom right panel,
where we plot the L-curve calculated from the normal equations.
The curve agrees with the TSQR curve for larger damping parameters,
but for small :math:`\lambda`, the normal equations approach cannot
provide accurate solution vectors leading to numerical
inaccuracies in the left portion of the curve.

.. _fig_multilarge:

.. figure:: /images/multilarge.png

   Top left: unregularized solutions; top right: regularized solutions;
   bottom left: L-curve for TSQR method; bottom right: L-curve from normal equations method.

.. include:: examples/largefit.c
   :code:

References and Further Reading
==============================

A summary of formulas and techniques for least squares fitting can be
found in the "Statistics" chapter of the Annual Review of Particle
Physics prepared by the Particle Data Group,

* *Review of Particle Properties*,
  R.M. Barnett et al., Physical Review D54, 1 (1996)
  http://pdg.lbl.gov

The Review of Particle Physics is available online at the website given
above.

.. index::
   single: NIST Statistical Reference Datasets
   single: Statistical Reference Datasets (StRD)

The tests used to prepare these routines are based on the NIST
Statistical Reference Datasets. The datasets and their documentation are
available from NIST at the following website,

http://www.nist.gov/itl/div898/strd/index.html

More information on Tikhonov regularization can be found in

* Hansen, P. C. (1998), Rank-Deficient and Discrete Ill-Posed Problems:
  Numerical Aspects of Linear Inversion. SIAM Monogr. on Mathematical
  Modeling and Computation, Society for Industrial and Applied Mathematics

* M. Rezghi and S. M. Hosseini (2009), A new variant of L-curve for
  Tikhonov regularization, Journal of Computational and Applied Mathematics,
  Volume 231, Issue 2, pages 914-924.

The GSL implementation of robust linear regression closely follows the publications

* DuMouchel, W. and F. O'Brien (1989), "Integrating a robust
  option into a multiple regression computing environment,"
  Computer Science and Statistics:  Proceedings of the 21st
  Symposium on the Interface, American Statistical Association

* Street, J.O., R.J. Carroll, and D. Ruppert (1988), "A note on
  computing robust regression estimates via iteratively
  reweighted least squares," The American Statistician, v. 42, 
  pp. 152-154.

More information about the normal equations and TSQR approach for solving
large linear least squares systems can be found in the publications

* Trefethen, L. N. and Bau, D. (1997), "Numerical Linear Algebra", SIAM.

* Demmel, J., Grigori, L., Hoemmen, M. F., and Langou, J.
  "Communication-optimal parallel and sequential QR and LU factorizations",
  UCB Technical Report No. UCB/EECS-2008-89, 2008.
