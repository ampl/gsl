.. index::
   single: random number distributions
   single: cumulative distribution functions (CDFs)
   single: CDFs, cumulative distribution functions
   single: inverse cumulative distribution functions
   single: quantile functions

.. _chap_random-number-distributions:

***************************
Random Number Distributions
***************************

.. include:: include.rst

This chapter describes functions for generating random variates and
computing their probability distributions.  Samples from the
distributions described in this chapter can be obtained using any of the
random number generators in the library as an underlying source of
randomness.  

In the simplest cases a non-uniform distribution can be obtained
analytically from the uniform distribution of a random number generator
by applying an appropriate transformation.  This method uses one call to
the random number generator.  More complicated distributions are created
by the *acceptance-rejection* method, which compares the desired
distribution against a distribution which is similar and known
analytically.  This usually requires several samples from the generator.

The library also provides cumulative distribution functions and inverse
cumulative distribution functions, sometimes referred to as quantile
functions.  The cumulative distribution functions and their inverses are
computed separately for the upper and lower tails of the distribution,
allowing full accuracy to be retained for small results.

The functions for random variates and probability density functions
described in this section are declared in :file:`gsl_randist.h`.  The
corresponding cumulative distribution functions are declared in
:file:`gsl_cdf.h`.

Note that the discrete random variate functions always
return a value of type :code:`unsigned int`, and on most platforms this
has a maximum value of

.. only:: not texinfo
   
   .. math:: 2^{32}-1 \approx 4.29 \times 10^9

.. only:: texinfo

   ::

      2^32-1 ~=~ 4.29e9

They should only be called with
a safe range of parameters (where there is a negligible probability of
a variate exceeding this limit) to prevent incorrect results due to
overflow.

Introduction
============

Continuous random number distributions are defined by a probability
density function, :math:`p(x)`, such that the probability of :math:`x`
occurring in the infinitesimal range :math:`x` to :math:`x + dx` is
:math:`p(x) dx`.

The cumulative distribution function for the lower tail :math:`P(x)` is
defined by the integral,

.. math:: P(x) = \int_{-\infty}^{x} dx' p(x')

and gives the probability of a variate taking a value less than :math:`x`.

The cumulative distribution function for the upper tail :math:`Q(x)` is
defined by the integral,

.. math:: Q(x) = \int_{x}^{+\infty} dx' p(x')

and gives the probability of a variate taking a value greater than :math:`x`.

The upper and lower cumulative distribution functions are related by
:math:`P(x) + Q(x) = 1` and satisfy :math:`0 \le P(x) \le 1`,
:math:`0 \le Q(x) \le 1`.

The inverse cumulative distributions, :math:`x = P^{-1}(P)`
and :math:`x = Q^{-1}(Q)`
give the values of :math:`x`
which correspond to a specific value of :math:`P` or :math:`Q`.  
They can be used to find confidence limits from probability values.

For discrete distributions the probability of sampling the integer
value :math:`k` is given by :math:`p(k)`, where :math:`\sum_k p(k) = 1`.
The cumulative distribution for the lower tail :math:`P(k)` of a
discrete distribution is defined as,

.. math:: P(k) = \sum_{i \le k} p(i) 

where the sum is over the allowed range of the distribution less than
or equal to :math:`k`.  

The cumulative distribution for the upper tail of a discrete
distribution :math:`Q(k)` is defined as

.. math:: Q(k) = \sum_{i > k} p(i) 

giving the sum of probabilities for all values greater than :math:`k`.
These two definitions satisfy the identity :math:`P(k)+Q(k)=1`.

If the range of the distribution is 1 to :math:`n` inclusive then
:math:`P(n) = 1`, :math:`Q(n) = 0` while :math:`P(1) = p(1)`,
:math:`Q(1) = 1 - p(1)`.

|newpage|

The Gaussian Distribution
=========================

.. index:: Gaussian distribution

.. function:: double gsl_ran_gaussian (const gsl_rng * r, double sigma)

   This function returns a Gaussian random variate, with mean zero and
   standard deviation :data:`sigma`.  The probability distribution for
   Gaussian random variates is,

   .. math:: p(x) dx = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2) dx

   for :math:`x` in the range :math:`-\infty` to :math:`+\infty`.  Use the
   transformation :math:`z = \mu + x` on the numbers returned by
   :func:`gsl_ran_gaussian` to obtain a Gaussian distribution with mean
   :math:`\mu`.  This function uses the Box-Muller algorithm which requires two
   calls to the random number generator :data:`r`.

.. function:: double gsl_ran_gaussian_pdf (double x, double sigma)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Gaussian distribution with standard deviation :data:`sigma`, using
   the formula given above.

   .. image:: /images/rand-gaussian.png

.. index:: Ziggurat method

.. function:: double gsl_ran_gaussian_ziggurat (const gsl_rng * r, double sigma)
              double gsl_ran_gaussian_ratio_method (const gsl_rng * r, double sigma)

   This function computes a Gaussian random variate using the alternative
   Marsaglia-Tsang ziggurat and Kinderman-Monahan-Leva ratio methods.  The
   Ziggurat algorithm is the fastest available algorithm in most cases.

.. function:: double gsl_ran_ugaussian (const gsl_rng * r)
              double gsl_ran_ugaussian_pdf (double x)
              double gsl_ran_ugaussian_ratio_method (const gsl_rng * r)

   These functions compute results for the unit Gaussian distribution.  They
   are equivalent to the functions above with a standard deviation of one,
   :data:`sigma` = 1.

.. function:: double gsl_cdf_gaussian_P (double x, double sigma)
              double gsl_cdf_gaussian_Q (double x, double sigma)
              double gsl_cdf_gaussian_Pinv (double P, double sigma)
              double gsl_cdf_gaussian_Qinv (double Q, double sigma)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the Gaussian
   distribution with standard deviation :data:`sigma`.

.. function:: double gsl_cdf_ugaussian_P (double x)
              double gsl_cdf_ugaussian_Q (double x)
              double gsl_cdf_ugaussian_Pinv (double P)
              double gsl_cdf_ugaussian_Qinv (double Q)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the unit Gaussian
   distribution.

|newpage|

The Gaussian Tail Distribution
==============================

.. index:: Gaussian Tail distribution

.. function:: double gsl_ran_gaussian_tail (const gsl_rng * r, double a, double sigma)

   This function provides random variates from the upper tail of a Gaussian
   distribution with standard deviation :data:`sigma`.  The values returned
   are larger than the lower limit :data:`a`, which must be positive.  The
   method is based on Marsaglia's famous rectangle-wedge-tail algorithm (Ann. 
   Math. Stat. 32, 894--899 (1961)), with this aspect explained in Knuth, v2,
   3rd ed, p139,586 (exercise 11).

   The probability distribution for Gaussian tail random variates is,

   .. math:: p(x) dx = {1 \over N(a;\sigma) \sqrt{2 \pi \sigma^2}} \exp (- x^2 / 2\sigma^2) dx

   for :math:`x > a` where :math:`N(a;\sigma)` is the normalization constant,

   .. only:: not texinfo

      .. math:: N(a;\sigma) = {1 \over 2} \hbox{erfc}\left({a \over \sqrt{2 \sigma^2}}\right).

   .. only:: texinfo

      ::

         N(a;\sigma) = (1/2) erfc(a / sqrt(2 sigma^2)).

.. function:: double gsl_ran_gaussian_tail_pdf (double x, double a, double sigma)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Gaussian tail distribution with standard deviation :data:`sigma` and
   lower limit :data:`a`, using the formula given above.

   .. image:: /images/rand-gaussian-tail.png

.. function:: double gsl_ran_ugaussian_tail (const gsl_rng * r, double a)
              double gsl_ran_ugaussian_tail_pdf (double x, double a)

   These functions compute results for the tail of a unit Gaussian
   distribution.  They are equivalent to the functions above with a standard
   deviation of one, :data:`sigma` = 1.

|newpage|

The Bivariate Gaussian Distribution
===================================

.. index::
   single: Bivariate Gaussian distribution
   single: two dimensional Gaussian distribution
   single: Gaussian distribution, bivariate

.. function:: void gsl_ran_bivariate_gaussian (const gsl_rng * r, double sigma_x, double sigma_y, double rho, double * x, double * y)

   This function generates a pair of correlated Gaussian variates, with
   mean zero, correlation coefficient :data:`rho` and standard deviations
   :data:`sigma_x` and :data:`sigma_y` in the :math:`x` and :math:`y` directions.
   The probability distribution for bivariate Gaussian random variates is,

   .. only:: not texinfo

      .. math:: p(x,y) dx dy = {1 \over 2 \pi \sigma_x \sigma_y \sqrt{1-\rho^2}} \exp \left(-{(x^2/\sigma_x^2 + y^2/\sigma_y^2 - 2 \rho x y/(\sigma_x\sigma_y)) \over 2(1-\rho^2)}\right) dx dy

   .. only:: texinfo

      ::

         p(x,y) dx dy = {1 \over 2 \pi \sigma_x \sigma_y \sqrt{1-\rho^2}} \exp (-(x^2/\sigma_x^2 + y^2/\sigma_y^2 - 2 \rho x y/(\sigma_x\sigma_y))/2(1-\rho^2)) dx dy

   for :math:`x,y` in the range :math:`-\infty` to :math:`+\infty`.  The
   correlation coefficient :data:`rho` should lie between :math:`1` and
   :math:`-1`.

.. function:: double gsl_ran_bivariate_gaussian_pdf (double x, double y, double sigma_x, double sigma_y, double rho)

   This function computes the probability density :math:`p(x,y)` at
   (:data:`x`, :data:`y`) for a bivariate Gaussian distribution with standard
   deviations :data:`sigma_x`, :data:`sigma_y` and correlation coefficient
   :data:`rho`, using the formula given above.

   .. image:: /images/rand-bivariate-gaussian.png

|newpage|

The Multivariate Gaussian Distribution
======================================

.. index::
   single: Bivariate Gaussian distribution
   single: two dimensional Gaussian distribution
   single: Gaussian distribution, bivariate

.. function:: int gsl_ran_multivariate_gaussian (const gsl_rng * r, const gsl_vector * mu, const gsl_matrix * L, gsl_vector * result)

   This function generates a random vector satisfying the :math:`k`-dimensional multivariate Gaussian
   distribution with mean :math:`\mu` and variance-covariance matrix
   :math:`\Sigma`. On input, the :math:`k`-vector :math:`\mu` is given in :data:`mu`, and
   the Cholesky factor of the :math:`k`-by-:math:`k` matrix :math:`\Sigma = L L^T` is
   given in the lower triangle of :data:`L`, as output from :func:`gsl_linalg_cholesky_decomp`.
   The random vector is stored in :data:`result` on output. The probability distribution
   for multivariate Gaussian random variates is

   .. only:: not texinfo

      .. math:: p(x_1,\dots,x_k) dx_1 \dots dx_k = {1 \over \sqrt{(2 \pi)^k |\Sigma|}} \exp \left(-{1 \over 2} (x - \mu)^T \Sigma^{-1} (x - \mu)\right) dx_1 \dots dx_k

   .. only:: texinfo

      ::

         p(x_1,...,x_k) dx_1 ... dx_k = 1 / ( \sqrt{(2 \pi)^k |\Sigma| ) \exp (-1/2 (x - \mu)^T \Sigma^{-1} (x - \mu)) dx_1 ... dx_k

.. function:: int gsl_ran_multivariate_gaussian_pdf (const gsl_vector * x, const gsl_vector * mu, const gsl_matrix * L, double * result, gsl_vector * work)
              int gsl_ran_multivariate_gaussian_log_pdf (const gsl_vector * x, const gsl_vector * mu, const gsl_matrix * L, double * result, gsl_vector * work)

   These functions compute :math:`p(x)` or :math:`\log{p(x)}` at the point :data:`x`, using mean vector
   :data:`mu` and variance-covariance matrix specified by its Cholesky factor :data:`L` using the formula
   above. Additional workspace of length :math:`k` is required in :data:`work`.

.. function:: int gsl_ran_multivariate_gaussian_mean (const gsl_matrix * X, gsl_vector * mu_hat)

   Given a set of :math:`n` samples :math:`X_j` from a :math:`k`-dimensional multivariate Gaussian distribution,
   this function computes the maximum likelihood estimate of the mean of the distribution, given by

   .. math:: \Hat{\mu} = {1 \over n} \sum_{j=1}^n X_j

   The samples :math:`X_1,X_2,\dots,X_n` are given in the :math:`n`-by-:math:`k` matrix :data:`X`, and the maximum
   likelihood estimate of the mean is stored in :data:`mu_hat` on output.

.. function:: int gsl_ran_multivariate_gaussian_vcov (const gsl_matrix * X, gsl_matrix * sigma_hat)

   Given a set of :math:`n` samples :math:`X_j` from a :math:`k`-dimensional multivariate Gaussian distribution,
   this function computes the maximum likelihood estimate of the variance-covariance matrix of the distribution,
   given by

   .. only:: not texinfo

      .. math:: \Hat{\Sigma} = {1 \over n} \sum_{j=1}^n \left( X_j - \Hat{\mu} \right) \left( X_j - \Hat{\mu} \right)^T

   .. only:: texinfo

      ::

         \Hat{\Sigma} = (1 / n) \sum_{j=1}^n ( X_j - \Hat{\mu} ) ( X_j - \Hat{\mu} )^T

   The samples :math:`X_1,X_2,\dots,X_n` are given in the :math:`n`-by-:math:`k` matrix :data:`X` and the maximum
   likelihood estimate of the variance-covariance matrix is stored in :data:`sigma_hat` on output.

|newpage|

The Exponential Distribution
============================

.. index:: Exponential distribution

.. function:: double gsl_ran_exponential (const gsl_rng * r, double mu)

   This function returns a random variate from the exponential distribution
   with mean :data:`mu`. The distribution is,

   .. math:: p(x) dx = {1 \over \mu} \exp(-x/\mu) dx

   for :math:`x \ge 0`.

.. function:: double gsl_ran_exponential_pdf (double x, double mu)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for an exponential distribution with mean :data:`mu`, using the formula
   given above.

   .. image:: /images/rand-exponential.png

.. function:: double gsl_cdf_exponential_P (double x, double mu)
              double gsl_cdf_exponential_Q (double x, double mu)
              double gsl_cdf_exponential_Pinv (double P, double mu)
              double gsl_cdf_exponential_Qinv (double Q, double mu)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the exponential
   distribution with mean :data:`mu`.

|newpage|

The Laplace Distribution
========================

.. index::
   single: two-sided exponential distribution
   single: Laplace distribution

.. function:: double gsl_ran_laplace (const gsl_rng * r, double a)

   This function returns a random variate from the Laplace distribution
   with width :data:`a`.  The distribution is,

   .. math:: p(x) dx = {1 \over 2 a}  \exp(-|x/a|) dx

   for :math:`-\infty < x < \infty`.

.. function:: double gsl_ran_laplace_pdf (double x, double a)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Laplace distribution with width :data:`a`, using the formula
   given above.

   .. image:: /images/rand-laplace.png

.. function:: double gsl_cdf_laplace_P (double x, double a)
              double gsl_cdf_laplace_Q (double x, double a)
              double gsl_cdf_laplace_Pinv (double P, double a)
              double gsl_cdf_laplace_Qinv (double Q, double a)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the Laplace
   distribution with width :data:`a`.

|newpage|

The Exponential Power Distribution
==================================

.. index:: Exponential power distribution

.. function:: double gsl_ran_exppow (const gsl_rng * r, double a, double b)

   This function returns a random variate from the exponential power distribution
   with scale parameter :data:`a` and exponent :data:`b`.  The distribution is,

   .. math:: p(x) dx = {1 \over 2 a \Gamma(1+1/b)} \exp(-|x/a|^b) dx

   for :math:`x \ge 0`.
   For :math:`b = 1` this reduces to the Laplace
   distribution.  For :math:`b = 2` it has the same form as a Gaussian
   distribution, but with :math:`a = \sqrt{2} \sigma`.

.. function:: double gsl_ran_exppow_pdf (double x, double a, double b)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for an exponential power distribution with scale parameter :data:`a`
   and exponent :data:`b`, using the formula given above.

   .. image:: /images/rand-exppow.png

.. function:: double gsl_cdf_exppow_P (double x, double a, double b)
              double gsl_cdf_exppow_Q (double x, double a, double b)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` for the exponential power distribution with
   parameters :data:`a` and :data:`b`.

|newpage|

The Cauchy Distribution
=======================

.. index:: Cauchy distribution

.. function:: double gsl_ran_cauchy (const gsl_rng * r, double a)

   This function returns a random variate from the Cauchy distribution with
   scale parameter :data:`a`.  The probability distribution for Cauchy
   random variates is,

   .. math:: p(x) dx = {1 \over a\pi (1 + (x/a)^2) } dx

   for :math:`x` in the range :math:`-\infty` to :math:`+\infty`.  The Cauchy
   distribution is also known as the Lorentz distribution.

.. function:: double gsl_ran_cauchy_pdf (double x, double a)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Cauchy distribution with scale parameter :data:`a`, using the formula
   given above.

   .. image:: /images/rand-cauchy.png

.. function:: double gsl_cdf_cauchy_P (double x, double a)
              double gsl_cdf_cauchy_Q (double x, double a)
              double gsl_cdf_cauchy_Pinv (double P, double a)
              double gsl_cdf_cauchy_Qinv (double Q, double a)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the Cauchy
   distribution with scale parameter :data:`a`.

|newpage|

The Rayleigh Distribution
=========================

.. index:: Rayleigh distribution

.. function:: double gsl_ran_rayleigh (const gsl_rng * r, double sigma)

   This function returns a random variate from the Rayleigh distribution with
   scale parameter :data:`sigma`.  The distribution is,

   .. math:: p(x) dx = {x \over \sigma^2} \exp(- x^2/(2 \sigma^2)) dx

   for :math:`x > 0`.

.. function:: double gsl_ran_rayleigh_pdf (double x, double sigma)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Rayleigh distribution with scale parameter :data:`sigma`, using the
   formula given above.

   .. image:: /images/rand-rayleigh.png

.. function:: double gsl_cdf_rayleigh_P (double x, double sigma)
              double gsl_cdf_rayleigh_Q (double x, double sigma)
              double gsl_cdf_rayleigh_Pinv (double P, double sigma)
              double gsl_cdf_rayleigh_Qinv (double Q, double sigma)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the Rayleigh
   distribution with scale parameter :data:`sigma`.

|newpage|

The Rayleigh Tail Distribution
==============================

.. index:: Rayleigh Tail distribution

.. function:: double gsl_ran_rayleigh_tail (const gsl_rng * r, double a, double sigma)

   This function returns a random variate from the tail of the Rayleigh
   distribution with scale parameter :data:`sigma` and a lower limit of
   :data:`a`.  The distribution is,

   .. math:: p(x) dx = {x \over \sigma^2} \exp ((a^2 - x^2) /(2 \sigma^2)) dx

   for :math:`x > a`.

.. function:: double gsl_ran_rayleigh_tail_pdf (double x, double a, double sigma)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Rayleigh tail distribution with scale parameter :data:`sigma` and
   lower limit :data:`a`, using the formula given above.

   .. image:: /images/rand-rayleigh-tail.png

|newpage|

The Landau Distribution
=======================

.. index:: Landau distribution

.. function:: double gsl_ran_landau (const gsl_rng * r)

   This function returns a random variate from the Landau distribution.  The
   probability distribution for Landau random variates is defined
   analytically by the complex integral,

   .. only:: not texinfo

      .. math:: p(x) = {1 \over {2 \pi i}} \int_{c-i\infty}^{c+i\infty} ds\, \exp(s \log(s) + x s) 

   .. only:: texinfo

      ::

         p(x) = (1/(2 \pi i)) \int_{c-i\infty}^{c+i\infty} ds exp(s log(s) + x s) 

   For numerical purposes it is more convenient to use the following
   equivalent form of the integral,

   .. math:: p(x) = (1/\pi) \int_0^\infty dt \exp(-t \log(t) - x t) \sin(\pi t).

.. function:: double gsl_ran_landau_pdf (double x)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for the Landau distribution using an approximation to the formula given
   above.

   .. image:: /images/rand-landau.png

|newpage|

The Levy alpha-Stable Distributions
===================================

.. index:: Levy distribution

.. function:: double gsl_ran_levy (const gsl_rng * r, double c, double alpha)

   This function returns a random variate from the Levy symmetric stable
   distribution with scale :data:`c` and exponent :data:`alpha`.  The symmetric
   stable probability distribution is defined by a Fourier transform,

   .. only:: not texinfo

      .. math:: p(x) = {1 \over 2 \pi} \int_{-\infty}^{+\infty} dt \exp(-it x - |c t|^\alpha)

   .. only:: texinfo

      ::

         p(x) = 1 / (2 \pi) \int_{-\infty}^{+\infty} dt \exp(-it x - |c t|^alpha)

   There is no explicit solution for the form of :math:`p(x)` and the
   library does not define a corresponding :code:`pdf` function.  For
   :math:`\alpha = 1` the distribution reduces to the Cauchy distribution.  For
   :math:`\alpha = 2` it is a Gaussian distribution with :math:`\sigma = \sqrt{2} c`.
   For :math:`\alpha < 1` the tails of the distribution become extremely wide.

   The algorithm only works for :math:`0 < \alpha \le 2`.

   .. image:: /images/rand-levy.png

|newpage|

The Levy skew alpha-Stable Distribution
=======================================

.. index::
   single: Levy distribution, skew
   single: Skew Levy distribution

.. function:: double gsl_ran_levy_skew (const gsl_rng * r, double c, double alpha, double beta)

   This function returns a random variate from the Levy skew stable
   distribution with scale :data:`c`, exponent :data:`alpha` and skewness
   parameter :data:`beta`.  The skewness parameter must lie in the range
   :math:`[-1,1]`.  The Levy skew stable probability distribution is defined
   by a Fourier transform,

   .. only:: not texinfo

      .. math:: p(x) = {1 \over 2 \pi} \int_{-\infty}^{+\infty} dt \exp(-it x - |c t|^\alpha (1-i \beta \sgn(t) \tan(\pi\alpha/2)))

   .. only:: texinfo

      ::

         p(x) = 1 / (2 \pi) \int_{-\infty}^{+\infty} dt \exp(-it x - |c t|^alpha (1-i beta sign(t) tan(pi alpha/2)))

   When :math:`\alpha = 1` the term :math:`\tan(\pi \alpha/2)` is replaced by
   :math:`-(2/\pi)\log|t|`.  There is no explicit solution for the form of
   :math:`p(x)` and the library does not define a corresponding :code:`pdf`
   function.  For :math:`\alpha = 2` the distribution reduces to a Gaussian
   distribution with :math:`\sigma = \sqrt{2} c`
   and the skewness parameter has no effect.  
   For :math:`\alpha < 1` the tails of the distribution become extremely
   wide.  The symmetric distribution corresponds to :math:`\beta = 0`.

   The algorithm only works for :math:`0 < \alpha \le 2`.

The Levy alpha-stable distributions have the property that if :math:`N`
alpha-stable variates are drawn from the distribution :math:`p(c, \alpha, \beta)`
then the sum :math:`Y = X_1 + X_2 + \dots + X_N` will also be
distributed as an alpha-stable variate,
:math:`p(N^{1/\alpha} c, \alpha, \beta)`.

.. PDF not available because there is no analytic expression for it
..
.. @deftypefun double gsl_ran_levy_pdf (double x, double mu)
.. This function computes the probability density :math:`p(x)` at :data:`x`
.. for a symmetric Levy distribution with scale parameter :data:`mu` and
.. exponent :data:`a`, using the formula given above.
.. @end deftypefun

.. image:: /images/rand-levyskew.png

|newpage|

The Gamma Distribution
======================

.. index::
   single: Gamma distribution
   single: Erlang distribution

.. function:: double gsl_ran_gamma (const gsl_rng * r, double a, double b)

   This function returns a random variate from the gamma
   distribution.  The distribution function is,

   .. math:: p(x) dx = {1 \over \Gamma(a) b^a} x^{a-1} e^{-x/b} dx

   for :math:`x > 0`.

   The gamma distribution with an integer parameter :data:`a` is known as the Erlang distribution.

   The variates are computed using the Marsaglia-Tsang fast gamma method.
   This function for this method was previously called
   :func:`gsl_ran_gamma_mt` and can still be accessed using this name.

.. If @xmath{X} and @xmath{Y} are independent gamma-distributed random
.. variables of order @xmath{a} and @xmath{b}, then @xmath{X+Y} has a gamma
.. distribution of order @xmath{a+b}.

.. function:: double gsl_ran_gamma_knuth (const gsl_rng * r, double a, double b)

   This function returns a gamma variate using the algorithms from Knuth (vol 2).

.. function:: double gsl_ran_gamma_pdf (double x, double a, double b)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a gamma distribution with parameters :data:`a` and :data:`b`, using the
   formula given above.

   .. image:: /images/rand-gamma.png

.. function:: double gsl_cdf_gamma_P (double x, double a, double b)
              double gsl_cdf_gamma_Q (double x, double a, double b)
              double gsl_cdf_gamma_Pinv (double P, double a, double b)
              double gsl_cdf_gamma_Qinv (double Q, double a, double b)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the gamma
   distribution with parameters :data:`a` and :data:`b`.

|newpage|

The Flat (Uniform) Distribution
===============================

.. index::
   single: flat distribution
   single: uniform distribution

.. function:: double gsl_ran_flat (const gsl_rng * r, double a, double b)

   This function returns a random variate from the flat (uniform)
   distribution from :data:`a` to :data:`b`. The distribution is,

   .. math:: p(x) dx = {1 \over (b-a)} dx

   if :math:`a \le x < b` and 0 otherwise.

.. function:: double gsl_ran_flat_pdf (double x, double a, double b)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a uniform distribution from :data:`a` to :data:`b`, using the formula
   given above.

   .. image:: /images/rand-flat.png

.. function:: double gsl_cdf_flat_P (double x, double a, double b)
              double gsl_cdf_flat_Q (double x, double a, double b)
              double gsl_cdf_flat_Pinv (double P, double a, double b)
              double gsl_cdf_flat_Qinv (double Q, double a, double b)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for a uniform distribution
   from :data:`a` to :data:`b`.

|newpage|

The Lognormal Distribution
==========================

.. index:: Lognormal distribution

.. function:: double gsl_ran_lognormal (const gsl_rng * r, double zeta, double sigma)

   This function returns a random variate from the lognormal
   distribution.  The distribution function is,

   .. math:: p(x) dx = {1 \over x \sqrt{2 \pi \sigma^2}} \exp(-(\ln(x) - \zeta)^2/2 \sigma^2) dx

   for :math:`x > 0`.

.. function:: double gsl_ran_lognormal_pdf (double x, double zeta, double sigma)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a lognormal distribution with parameters :data:`zeta` and :data:`sigma`,
   using the formula given above.

   .. image:: /images/rand-lognormal.png

.. function:: double gsl_cdf_lognormal_P (double x, double zeta, double sigma)
              double gsl_cdf_lognormal_Q (double x, double zeta, double sigma)
              double gsl_cdf_lognormal_Pinv (double P, double zeta, double sigma)
              double gsl_cdf_lognormal_Qinv (double Q, double zeta, double sigma)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the lognormal
   distribution with parameters :data:`zeta` and :data:`sigma`.

|newpage|

The Chi-squared Distribution
============================

The chi-squared distribution arises in statistics.  If :math:`Y_i` are
:math:`n` independent Gaussian random variates with unit variance then the
sum-of-squares,

.. math:: X_i = \sum_i Y_i^2

has a chi-squared distribution with :math:`n` degrees of freedom.

.. index:: Chi-squared distribution

.. function:: double gsl_ran_chisq (const gsl_rng * r, double nu)

   This function returns a random variate from the chi-squared distribution
   with :data:`nu` degrees of freedom. The distribution function is,

   .. math:: p(x) dx = {1 \over 2 \Gamma(\nu/2) } (x/2)^{\nu/2 - 1} \exp(-x/2) dx

   for :math:`x \ge 0`.

.. function:: double gsl_ran_chisq_pdf (double x, double nu)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a chi-squared distribution with :data:`nu` degrees of freedom, using
   the formula given above.

   .. image:: /images/rand-chisq.png

.. function:: double gsl_cdf_chisq_P (double x, double nu)
              double gsl_cdf_chisq_Q (double x, double nu)
              double gsl_cdf_chisq_Pinv (double P, double nu)
              double gsl_cdf_chisq_Qinv (double Q, double nu)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the chi-squared
   distribution with :data:`nu` degrees of freedom.

|newpage|

The F-distribution
==================

The F-distribution arises in statistics.  If :math:`Y_1` and :math:`Y_2`
are chi-squared deviates with :math:`\nu_1` and :math:`\nu_2` degrees of
freedom then the ratio,

.. math:: X = { (Y_1 / \nu_1) \over (Y_2 / \nu_2) }

has an F-distribution :math:`F(x;\nu_1,\nu_2)`.

.. index:: F-distribution

.. function:: double gsl_ran_fdist (const gsl_rng * r, double nu1, double nu2)

   This function returns a random variate from the F-distribution with degrees of freedom :data:`nu1` and :data:`nu2`.
   The distribution function is,

   .. only:: not texinfo

      .. math::

         p(x) dx = 
            { \Gamma((\nu_1 + \nu_2)/2)
                 \over \Gamma(\nu_1/2) \Gamma(\nu_2/2) } 
            \nu_1^{\nu_1/2} \nu_2^{\nu_2/2} 
            x^{\nu_1/2 - 1} (\nu_2 + \nu_1 x)^{-\nu_1/2 -\nu_2/2}

   .. only:: texinfo

      ::

         p(x) dx = 
            { \Gamma((\nu_1 + \nu_2)/2)
                 \over \Gamma(\nu_1/2) \Gamma(\nu_2/2) } 
            \nu_1^{\nu_1/2} \nu_2^{\nu_2/2} 
            x^{\nu_1/2 - 1} (\nu_2 + \nu_1 x)^{-\nu_1/2 -\nu_2/2}

   for :math:`x \ge 0`.

.. function:: double gsl_ran_fdist_pdf (double x, double nu1, double nu2)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for an F-distribution with :data:`nu1` and :data:`nu2` degrees of freedom,
   using the formula given above.

   .. image:: /images/rand-fdist.png

.. function:: double gsl_cdf_fdist_P (double x, double nu1, double nu2)
              double gsl_cdf_fdist_Q (double x, double nu1, double nu2)
              double gsl_cdf_fdist_Pinv (double P, double nu1, double nu2)
              double gsl_cdf_fdist_Qinv (double Q, double nu1, double nu2)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the F-distribution
   with :data:`nu1` and :data:`nu2` degrees of freedom.

|newpage|

The t-distribution
==================

The t-distribution arises in statistics.  If :math:`Y_1` has a normal
distribution and :math:`Y_2` has a chi-squared distribution with
:math:`\nu` degrees of freedom then the ratio,

.. math:: X = { Y_1 \over \sqrt{Y_2 / \nu} }

has a t-distribution :math:`t(x;\nu)` with :math:`\nu` degrees of freedom.

.. index::
   single: t-distribution
   single: Student t-distribution

.. function:: double gsl_ran_tdist (const gsl_rng * r, double nu)

   This function returns a random variate from the t-distribution.  The
   distribution function is,

   .. math::

      p(x) dx = {\Gamma((\nu + 1)/2) \over \sqrt{\pi \nu} \Gamma(\nu/2)}
         (1 + x^2/\nu)^{-(\nu + 1)/2} dx

   for :math:`-\infty < x < +\infty`.

.. function:: double gsl_ran_tdist_pdf (double x, double nu)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a t-distribution with :data:`nu` degrees of freedom, using the formula
   given above.

   .. image:: /images/rand-tdist.png

.. function:: double gsl_cdf_tdist_P (double x, double nu)
              double gsl_cdf_tdist_Q (double x, double nu)
              double gsl_cdf_tdist_Pinv (double P, double nu)
              double gsl_cdf_tdist_Qinv (double Q, double nu)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the t-distribution
   with :data:`nu` degrees of freedom.

|newpage|

The Beta Distribution
=====================

.. index:: Beta distribution

.. function:: double gsl_ran_beta (const gsl_rng * r, double a, double b)

   This function returns a random variate from the beta
   distribution.  The distribution function is,

   .. math:: p(x) dx = {\Gamma(a+b) \over \Gamma(a) \Gamma(b)} x^{a-1} (1-x)^{b-1} dx

   for :math:`0 \le x \le 1`.

.. function:: double gsl_ran_beta_pdf (double x, double a, double b)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a beta distribution with parameters :data:`a` and :data:`b`, using the
   formula given above.

   .. image:: /images/rand-beta.png

.. function:: double gsl_cdf_beta_P (double x, double a, double b)
              double gsl_cdf_beta_Q (double x, double a, double b)
              double gsl_cdf_beta_Pinv (double P, double a, double b)
              double gsl_cdf_beta_Qinv (double Q, double a, double b)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the beta
   distribution with parameters :data:`a` and :data:`b`.

|newpage|

The Logistic Distribution
=========================

.. index:: Logistic distribution

.. function:: double gsl_ran_logistic (const gsl_rng * r, double a)

   This function returns a random variate from the logistic
   distribution.  The distribution function is,

   .. math:: p(x) dx = { \exp(-x/a) \over a (1 + \exp(-x/a))^2 } dx

   for :math:`-\infty < x < +\infty`.

.. function:: double gsl_ran_logistic_pdf (double x, double a)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a logistic distribution with scale parameter :data:`a`, using the
   formula given above.

   .. image:: /images/rand-logistic.png

.. function:: double gsl_cdf_logistic_P (double x, double a)
              double gsl_cdf_logistic_Q (double x, double a)
              double gsl_cdf_logistic_Pinv (double P, double a)
              double gsl_cdf_logistic_Qinv (double Q, double a)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the logistic
   distribution with scale parameter :data:`a`.

|newpage|

The Pareto Distribution
=======================

.. index:: Pareto distribution

.. function:: double gsl_ran_pareto (const gsl_rng * r, double a, double b)

   This function returns a random variate from the Pareto distribution of
   order :data:`a`.  The distribution function is,

   .. math:: p(x) dx = (a/b) / (x/b)^{a+1} dx

   for :math:`x \ge b`.

.. function:: double gsl_ran_pareto_pdf (double x, double a, double b)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Pareto distribution with exponent :data:`a` and scale :data:`b`, using
   the formula given above.

   .. image:: /images/rand-pareto.png

.. function:: double gsl_cdf_pareto_P (double x, double a, double b)
              double gsl_cdf_pareto_Q (double x, double a, double b)
              double gsl_cdf_pareto_Pinv (double P, double a, double b)
              double gsl_cdf_pareto_Qinv (double Q, double a, double b)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the Pareto
   distribution with exponent :data:`a` and scale :data:`b`.

|newpage|

Spherical Vector Distributions
==============================

The spherical distributions generate random vectors, located on a
spherical surface.  They can be used as random directions, for example in
the steps of a random walk.

.. index::
   single: 2D random direction vector
   single: direction vector, random 2D
   single: spherical random variates, 2D

.. function:: void gsl_ran_dir_2d (const gsl_rng * r, double * x, double * y)
              void gsl_ran_dir_2d_trig_method (const gsl_rng * r, double * x, double * y)

   This function returns a random direction vector :math:`v` =
   (:data:`x`, :data:`y`) in two dimensions.  The vector is normalized such that
   :math:`|v|^2 = x^2 + y^2 = 1`.  The obvious way to do this is to take a
   uniform random number between 0 and :math:`2\pi` and let :data:`x` and
   :data:`y` be the sine and cosine respectively.  Two trig functions would
   have been expensive in the old days, but with modern hardware
   implementations, this is sometimes the fastest way to go.  This is the
   case for the Pentium (but not the case for the Sun Sparcstation).
   One can avoid the trig evaluations by choosing :data:`x` and
   :data:`y` in the interior of a unit circle (choose them at random from the
   interior of the enclosing square, and then reject those that are outside
   the unit circle), and then dividing by :math:`\sqrt{x^2 + y^2}`.
   A much cleverer approach, attributed to von Neumann (See Knuth, v2, 3rd
   ed, p140, exercise 23), requires neither trig nor a square root.  In
   this approach, :data:`u` and :data:`v` are chosen at random from the
   interior of a unit circle, and then :math:`x=(u^2-v^2)/(u^2+v^2)` and
   :math:`y=2uv/(u^2+v^2)`.

.. index::
   single: 3D random direction vector
   single: direction vector, random 3D
   single: spherical random variates, 3D

.. function:: void gsl_ran_dir_3d (const gsl_rng * r, double * x, double * y, double * z)

   This function returns a random direction vector :math:`v` =
   (:data:`x`, :data:`y`, :data:`z`) in three dimensions.  The vector is normalized
   such that :math:`|v|^2 = x^2 + y^2 + z^2 = 1`.  The method employed is
   due to Robert E. Knop (CACM 13, 326 (1970)), and explained in Knuth, v2,
   3rd ed, p136.  It uses the surprising fact that the distribution
   projected along any axis is actually uniform (this is only true for 3
   dimensions).

.. index::
   single: N-dimensional random direction vector
   single: direction vector, random N-dimensional
   single: spherical random variates, N-dimensional

.. function:: void gsl_ran_dir_nd (const gsl_rng * r, size_t n, double * x)

   This function returns a random direction vector
   :math:`v = (x_1,x_2,\ldots,x_n)`
   in :data:`n` dimensions.  The vector is normalized such that 
   :math:`|v|^2 = x_1^2 + x_2^2 + \cdots + x_n^2 = 1`.
   The method
   uses the fact that a multivariate Gaussian distribution is spherically
   symmetric.  Each component is generated to have a Gaussian distribution,
   and then the components are normalized.  The method is described by
   Knuth, v2, 3rd ed, p135--136, and attributed to G. W. Brown, Modern
   Mathematics for the Engineer (1956).

|newpage|

The Weibull Distribution
========================

.. index:: Weibull distribution

.. function:: double gsl_ran_weibull (const gsl_rng * r, double a, double b)

   This function returns a random variate from the Weibull distribution.  The
   distribution function is,

   .. math:: p(x) dx = {b \over a^b} x^{b-1}  \exp(-(x/a)^b) dx

   for :math:`x \ge 0`.

.. function:: double gsl_ran_weibull_pdf (double x, double a, double b)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Weibull distribution with scale :data:`a` and exponent :data:`b`,
   using the formula given above.

   .. image:: /images/rand-weibull.png

.. function:: double gsl_cdf_weibull_P (double x, double a, double b)
              double gsl_cdf_weibull_Q (double x, double a, double b)
              double gsl_cdf_weibull_Pinv (double P, double a, double b)
              double gsl_cdf_weibull_Qinv (double Q, double a, double b)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the Weibull
   distribution with scale :data:`a` and exponent :data:`b`.

|newpage|

The Type-1 Gumbel Distribution
==============================

.. index::
   single: Gumbel distribution (Type 1)
   single: Type 1 Gumbel distribution, random variates

.. function:: double gsl_ran_gumbel1 (const gsl_rng * r, double a, double b)

   This function returns  a random variate from the Type-1 Gumbel
   distribution.  The Type-1 Gumbel distribution function is,

   .. math:: p(x) dx = a b \exp(-(b \exp(-ax) + ax)) dx

   for :math:`-\infty < x < \infty`. 

.. function:: double gsl_ran_gumbel1_pdf (double x, double a, double b)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Type-1 Gumbel distribution with parameters :data:`a` and :data:`b`,
   using the formula given above.

   .. image:: /images/rand-gumbel1.png

.. function:: double gsl_cdf_gumbel1_P (double x, double a, double b)
              double gsl_cdf_gumbel1_Q (double x, double a, double b)
              double gsl_cdf_gumbel1_Pinv (double P, double a, double b)
              double gsl_cdf_gumbel1_Qinv (double Q, double a, double b)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the Type-1 Gumbel
   distribution with parameters :data:`a` and :data:`b`.

|newpage|

The Type-2 Gumbel Distribution
==============================

.. index::
   single: Gumbel distribution (Type 2)
   single: Type 2 Gumbel distribution

.. function:: double gsl_ran_gumbel2 (const gsl_rng * r, double a, double b)

   This function returns a random variate from the Type-2 Gumbel
   distribution.  The Type-2 Gumbel distribution function is,

   .. math:: p(x) dx = a b x^{-a-1} \exp(-b x^{-a}) dx

   for :math:`0 < x < \infty`.

.. function:: double gsl_ran_gumbel2_pdf (double x, double a, double b)

   This function computes the probability density :math:`p(x)` at :data:`x`
   for a Type-2 Gumbel distribution with parameters :data:`a` and :data:`b`,
   using the formula given above.

   .. image:: /images/rand-gumbel2.png

.. function:: double gsl_cdf_gumbel2_P (double x, double a, double b)
              double gsl_cdf_gumbel2_Q (double x, double a, double b)
              double gsl_cdf_gumbel2_Pinv (double P, double a, double b)
              double gsl_cdf_gumbel2_Qinv (double Q, double a, double b)

   These functions compute the cumulative distribution functions
   :math:`P(x)`, :math:`Q(x)` and their inverses for the Type-2 Gumbel
   distribution with parameters :data:`a` and :data:`b`.

|newpage|

The Dirichlet Distribution
==========================

.. index:: Dirichlet distribution 

.. function:: void gsl_ran_dirichlet (const gsl_rng * r, size_t K, const double alpha[], double theta[])

   This function returns an array of :data:`K` random variates from a Dirichlet
   distribution of order :data:`K`-1. The distribution function is

   .. only:: not texinfo

      .. math::

         p(\theta_1,\ldots,\theta_K) \, d\theta_1 \cdots d\theta_K = 
                 {1 \over Z} \prod_{i=1}^{K} \theta_i^{\alpha_i - 1} 
                   \; \delta(1 -\sum_{i=1}^K \theta_i) d\theta_1 \cdots d\theta_K

   .. only:: texinfo

      ::

         p(\theta_1, ..., \theta_K) d\theta_1 ... d\theta_K = 
           (1/Z) \prod_{i=1}^K \theta_i^{\alpha_i - 1} \delta(1 -\sum_{i=1}^K \theta_i) d\theta_1 ... d\theta_K

   for :math:`\theta_i \ge 0`
   and :math:`\alpha_i > 0`.
   The delta function ensures that :math:`\sum \theta_i = 1`.
   The normalization factor :math:`Z` is

   .. math:: Z = {\prod_{i=1}^K \Gamma(\alpha_i) \over \Gamma( \sum_{i=1}^K \alpha_i)}

   The random variates are generated by sampling :data:`K` values 
   from gamma distributions with parameters 
   :math:`a=\alpha_i$, $b=1`,
   and renormalizing. 
   See A.M. Law, W.D. Kelton, *Simulation Modeling and Analysis* (1991).

.. function:: double gsl_ran_dirichlet_pdf (size_t K, const double alpha[], const double theta[]) 

   This function computes the probability density 
   :math:`p(\theta_1, \ldots , \theta_K)`
   at :code:`theta[K]` for a Dirichlet distribution with parameters 
   :code:`alpha[K]`, using the formula given above.

.. function:: double gsl_ran_dirichlet_lnpdf (size_t K, const double alpha[], const double theta[]) 

   This function computes the logarithm of the probability density 
   :math:`p(\theta_1, \ldots , \theta_K)`
   for a Dirichlet distribution with parameters 
   :code:`alpha[K]`.

|newpage|

General Discrete Distributions
==============================

Given :math:`K` discrete events with different probabilities :math:`P[k]`,
produce a random value :math:`k` consistent with its probability.

The obvious way to do this is to preprocess the probability list by
generating a cumulative probability array with :math:`K + 1` elements:

.. only:: not texinfo

   .. math::

      C[0] & = 0 \\
      C[k+1] &= C[k] + P[k]

.. only:: texinfo

   ::

     C[0] = 0 
     C[k+1] = C[k] + P[k]

Note that this construction produces :math:`C[K] = 1`.  Now choose a
uniform deviate :math:`u` between 0 and 1, and find the value of :math:`k`
such that :math:`C[k] \le u < C[k+1]`.
Although this in principle requires of order :math:`\log K` steps per
random number generation, they are fast steps, and if you use something
like :math:`\lfloor uK \rfloor` as a starting point, you can often do
pretty well.

But faster methods have been devised.  Again, the idea is to preprocess
the probability list, and save the result in some form of lookup table;
then the individual calls for a random discrete event can go rapidly.
An approach invented by G. Marsaglia (Generating discrete random variables
in a computer, Comm ACM 6, 37--38 (1963)) is very clever, and readers
interested in examples of good algorithm design are directed to this
short and well-written paper.  Unfortunately, for large :math:`K`,
Marsaglia's lookup table can be quite large.  

A much better approach is due to Alastair J. Walker (An efficient method
for generating discrete random variables with general distributions, ACM
Trans on Mathematical Software 3, 253--256 (1977); see also Knuth, v2,
3rd ed, p120--121,139).  This requires two lookup tables, one floating
point and one integer, but both only of size :math:`K`.  After
preprocessing, the random numbers are generated in O(1) time, even for
large :math:`K`.  The preprocessing suggested by Walker requires
:math:`O(K^2)` effort, but that is not actually necessary, and the
implementation provided here only takes :math:`O(K)` effort.  In general,
more preprocessing leads to faster generation of the individual random
numbers, but a diminishing return is reached pretty early.  Knuth points
out that the optimal preprocessing is combinatorially difficult for
large :math:`K`.

This method can be used to speed up some of the discrete random number
generators below, such as the binomial distribution.  To use it for
something like the Poisson Distribution, a modification would have to
be made, since it only takes a finite set of :math:`K` outcomes.

.. index::
   single: Discrete random numbers
   single: Discrete random numbers, preprocessing

.. type:: gsl_ran_discrete_t

   This structure contains the lookup table for the discrete random number generator.

.. function:: gsl_ran_discrete_t * gsl_ran_discrete_preproc (size_t K, const double * P)

   This function returns a pointer to a structure that contains the lookup
   table for the discrete random number generator.  The array :data:`P` contains
   the probabilities of the discrete events; these array elements must all be 
   positive, but they needn't add up to one (so you can think of them more
   generally as "weights")---the preprocessor will normalize appropriately.
   This return value is used
   as an argument for the :func:`gsl_ran_discrete` function below.

.. index:: Discrete random numbers

.. function:: size_t gsl_ran_discrete (const gsl_rng * r, const gsl_ran_discrete_t * g)

   After the preprocessor, above, has been called, you use this function to
   get the discrete random numbers.

.. index:: Discrete random numbers

.. function:: double gsl_ran_discrete_pdf (size_t k, const gsl_ran_discrete_t * g)

   Returns the probability :math:`P[k]` of observing the variable :data:`k`.
   Since :math:`P[k]` is not stored as part of the lookup table, it must be
   recomputed; this computation takes :math:`O(K)`, so if :data:`K` is large
   and you care about the original array :math:`P[k]` used to create the
   lookup table, then you should just keep this original array :math:`P[k]`
   around.

.. index:: Discrete random numbers

.. function:: void gsl_ran_discrete_free (gsl_ran_discrete_t * g)

   De-allocates the lookup table pointed to by :data:`g`.

|newpage|

The Poisson Distribution
========================

.. index:: Poisson random numbers

.. function:: unsigned int gsl_ran_poisson (const gsl_rng * r, double mu)

   This function returns a random integer from the Poisson distribution
   with mean :data:`mu`.  The probability distribution for Poisson variates is,

   .. math:: p(k) = {\mu^k \over k!} \exp(-\mu)

   for :math:`k \ge 0`.

.. function:: double gsl_ran_poisson_pdf (unsigned int k, double mu)

   This function computes the probability :math:`p(k)` of obtaining  :data:`k`
   from a Poisson distribution with mean :data:`mu`, using the formula
   given above.

   .. image:: /images/rand-poisson.png

.. function:: double gsl_cdf_poisson_P (unsigned int k, double mu)
              double gsl_cdf_poisson_Q (unsigned int k, double mu)

   These functions compute the cumulative distribution functions
   :math:`P(k)`, :math:`Q(k)` for the Poisson distribution with parameter
   :data:`mu`.

|newpage|

The Bernoulli Distribution
==========================

.. index::
   single: Bernoulli trial, random variates

.. function:: unsigned int gsl_ran_bernoulli (const gsl_rng * r, double p)

   This function returns either 0 or 1, the result of a Bernoulli trial
   with probability :data:`p`.  The probability distribution for a Bernoulli
   trial is,

   .. only:: not texinfo

      .. math::

         p(0) & = 1 - p \\
         p(1) & = p

   .. only:: texinfo

      ::

         p(0) = 1 - p
         p(1) = p

.. function:: double gsl_ran_bernoulli_pdf (unsigned int k, double p)

   This function computes the probability :math:`p(k)` of obtaining
   :data:`k` from a Bernoulli distribution with probability parameter
   :data:`p`, using the formula given above.

   .. image:: /images/rand-bernoulli.png

|newpage|

The Binomial Distribution
=========================

.. index:: Binomial random variates

.. function:: unsigned int gsl_ran_binomial (const gsl_rng * r, double p, unsigned int n)

   This function returns a random integer from the binomial distribution,
   the number of successes in :data:`n` independent trials with probability
   :data:`p`.  The probability distribution for binomial variates is,

   .. math:: p(k) = {n! \over k! (n-k)!} p^k (1-p)^{n-k}

   for :math:`0 \le k \le n`.

.. function:: double gsl_ran_binomial_pdf (unsigned int k, double p, unsigned int n)

   This function computes the probability :math:`p(k)` of obtaining :data:`k`
   from a binomial distribution with parameters :data:`p` and :data:`n`, using
   the formula given above.

   .. image:: /images/rand-binomial.png

.. function:: double gsl_cdf_binomial_P (unsigned int k, double p, unsigned int n)
              double gsl_cdf_binomial_Q (unsigned int k, double p, unsigned int n)

   These functions compute the cumulative distribution functions
   :math:`P(k)`, :math:`Q(k)`  for the binomial
   distribution with parameters :data:`p` and :data:`n`.

|newpage|

The Multinomial Distribution
============================

.. index:: Multinomial distribution 

.. function:: void gsl_ran_multinomial (const gsl_rng * r, size_t K, unsigned int N, const double p[], unsigned int n[])

   This function computes a random sample :data:`n` from the multinomial
   distribution formed by :data:`N` trials from an underlying distribution
   :code:`p[K]`. The distribution function for :data:`n` is,

   .. only:: not texinfo

      .. math::

         P(n_1, n_2,\cdots, n_K) = {{ N!}\over{n_1 ! n_2 ! \cdots n_K !}} \,
           p_1^{n_1} p_2^{n_2} \cdots p_K^{n_K}

   .. only:: texinfo

      ::

         P(n_1, n_2, ..., n_K) = 
           (N!/(n_1! n_2! ... n_K!)) p_1^n_1 p_2^n_2 ... p_K^n_K

   where :math:`(n_1, n_2, \ldots, n_K)`
   are nonnegative integers with 
   :math:`\sum_{k=1}^{K} n_k = N`,
   and
   :math:`(p_1, p_2, \ldots, p_K)`
   is a probability distribution with :math:`\sum p_i = 1`.  
   If the array :code:`p[K]` is not normalized then its entries will be
   treated as weights and normalized appropriately.  The arrays :data:`n`
   and :data:`p` must both be of length :data:`K`.

   Random variates are generated using the conditional binomial method (see
   C.S. Davis, *The computer generation of multinomial random
   variates*, Comp. Stat. Data Anal. 16 (1993) 205--217 for details).

.. function:: double gsl_ran_multinomial_pdf (size_t K, const double p[], const unsigned int n[]) 

   This function computes the probability 
   :math:`P(n_1, n_2, \ldots, n_K)`
   of sampling :code:`n[K]` from a multinomial distribution 
   with parameters :code:`p[K]`, using the formula given above.

.. function:: double gsl_ran_multinomial_lnpdf (size_t K, const double p[], const unsigned int n[]) 

   This function returns the logarithm of the probability for the
   multinomial distribution :math:`P(n_1, n_2, \ldots, n_K)`
   with parameters :code:`p[K]`.

|newpage|

The Negative Binomial Distribution
==================================

.. index::
   single: Negative Binomial distribution, random variates

.. function:: unsigned int gsl_ran_negative_binomial (const gsl_rng * r, double p, double n)

   This function returns a random integer from the negative binomial
   distribution, the number of failures occurring before :data:`n` successes
   in independent trials with probability :data:`p` of success.  The
   probability distribution for negative binomial variates is,

   .. math:: p(k) = {\Gamma(n + k) \over \Gamma(k+1) \Gamma(n) } p^n (1-p)^k

   Note that :math:`n` is not required to be an integer.

.. function:: double gsl_ran_negative_binomial_pdf (unsigned int k, double p, double  n)

   This function computes the probability :math:`p(k)` of obtaining :data:`k`
   from a negative binomial distribution with parameters :data:`p` and
   :data:`n`, using the formula given above.

   .. image:: /images/rand-nbinomial.png

.. function:: double gsl_cdf_negative_binomial_P (unsigned int k, double p, double n)
              double gsl_cdf_negative_binomial_Q (unsigned int k, double p, double n)

   These functions compute the cumulative distribution functions
   :math:`P(k)`, :math:`Q(k)` for the negative binomial distribution with
   parameters :data:`p` and :data:`n`.

|newpage|

The Pascal Distribution
=======================

.. function:: unsigned int gsl_ran_pascal (const gsl_rng * r, double p, unsigned int n)

   This function returns a random integer from the Pascal distribution.  The
   Pascal distribution is simply a negative binomial distribution with an
   integer value of :math:`n`.

   .. math:: p(k) = {(n + k - 1)! \over k! (n - 1)! } p^n (1-p)^k

   for :math:`k \ge 0`.

.. function:: double gsl_ran_pascal_pdf (unsigned int k, double p, unsigned int n)

   This function computes the probability :math:`p(k)` of obtaining :data:`k`
   from a Pascal distribution with parameters :data:`p` and
   :data:`n`, using the formula given above.

   .. image:: /images/rand-pascal.png

.. function:: double gsl_cdf_pascal_P (unsigned int k, double p, unsigned int n)
              double gsl_cdf_pascal_Q (unsigned int k, double p, unsigned int n)

   These functions compute the cumulative distribution functions
   :math:`P(k)`, :math:`Q(k)` for the Pascal distribution with
   parameters :data:`p` and :data:`n`.

|newpage|

The Geometric Distribution
==========================

.. index:: Geometric random variates

.. function:: unsigned int gsl_ran_geometric (const gsl_rng * r, double p)

   This function returns a random integer from the geometric distribution,
   the number of independent trials with probability :data:`p` until the
   first success.  The probability distribution for geometric variates
   is,

   .. math:: p(k) = p (1-p)^{k-1}

   for :math:`k \ge 1`.
   Note that the distribution begins with :math:`k = 1` with this
   definition.  There is another convention in which the exponent :math:`k - 1` 
   is replaced by :math:`k`.

.. function:: double gsl_ran_geometric_pdf (unsigned int k, double p)

   This function computes the probability :math:`p(k)` of obtaining :data:`k`
   from a geometric distribution with probability parameter :data:`p`, using
   the formula given above.

   .. image:: /images/rand-geometric.png

.. function:: double gsl_cdf_geometric_P (unsigned int k, double p)
              double gsl_cdf_geometric_Q (unsigned int k, double p)

   These functions compute the cumulative distribution functions
   :math:`P(k)`, :math:`Q(k)` for the geometric distribution with parameter
   :data:`p`.

|newpage|

The Hypergeometric Distribution
===============================

.. index::
   single: hypergeometric random variates
   single: Geometric random variates

.. function:: unsigned int gsl_ran_hypergeometric (const gsl_rng * r, unsigned int n1, unsigned int n2, unsigned int t)

   This function returns a random integer from the hypergeometric
   distribution.  The probability distribution for hypergeometric
   random variates is,

   .. math:: p(k) = C(n_1, k) C(n_2, t - k) / C(n_1 + n_2, t)

   where :math:`C(a,b) = a!/(b!(a-b)!)` and 
   :math:`t \leq n_1 + n_2`.
   The domain of :math:`k` is 
   :math:`\max(0, t - n_2), \ldots, \min(t, n_1)`

   If a population contains :math:`n_1` elements of "type 1" and
   :math:`n_2` elements of "type 2" then the hypergeometric
   distribution gives the probability of obtaining :math:`k` elements of
   "type 1" in :math:`t` samples from the population without
   replacement.

.. function:: double gsl_ran_hypergeometric_pdf (unsigned int k, unsigned int n1, unsigned int n2, unsigned int t)

   This function computes the probability :math:`p(k)` of obtaining :data:`k`
   from a hypergeometric distribution with parameters :data:`n1`, :data:`n2`,
   :data:`t`, using the formula given above.

   .. image:: /images/rand-hypergeometric.png

.. function:: double gsl_cdf_hypergeometric_P (unsigned int k, unsigned int n1, unsigned int n2, unsigned int t)
              double gsl_cdf_hypergeometric_Q (unsigned int k, unsigned int n1, unsigned int n2, unsigned int t)

   These functions compute the cumulative distribution functions
   :math:`P(k)`, :math:`Q(k)` for the hypergeometric distribution with
   parameters :data:`n1`, :data:`n2` and :data:`t`.

|newpage|

.. index:: Logarithmic random variates

The Logarithmic Distribution
============================

.. function:: unsigned int gsl_ran_logarithmic (const gsl_rng * r, double p)

   This function returns a random integer from the logarithmic
   distribution.  The probability distribution for logarithmic random variates
   is,

   .. only:: not texinfo

      .. math:: p(k) = {-1 \over \log(1-p)} {\left( p^k \over k \right)}

   .. only:: texinfo
      
      ::

         p(k) = {-1 \over \log(1-p)} {(p^k \over k)}

   for :math:`k \ge 1`.

.. function:: double gsl_ran_logarithmic_pdf (unsigned int k, double p)

   This function computes the probability :math:`p(k)` of obtaining :data:`k`
   from a logarithmic distribution with probability parameter :data:`p`,
   using the formula given above.

   .. image:: /images/rand-logarithmic.png

|newpage|

.. index:: Wishart random variates

The Wishart Distribution
========================

.. function:: int gsl_ran_wishart (const gsl_rng * r, const double n, const gsl_matrix * L, gsl_matrix * result, gsl_matrix * work)

   This function computes a random symmetric :math:`p`-by-:math:`p` matrix from the Wishart distribution.
   The probability distribution for Wishart random variates is,

   .. only:: not texinfo

      .. math:: p(X) = \frac{|X|^{(n-p-1)/2} e^{-\textrm{tr}\left( V^{-1} X\right)/2}}{2^{\frac{np}{2}} \left| V \right|^{n/2} \Gamma_p(\frac{n}{2})}

   .. only:: texinfo

      .. math:: p(X) = \frac{|X|^{(n-p-1)/2} e^{-tr( V^{-1} X)/2}}{2^{(np)/2} |V|^{n/2} \Gamma_p(n/2)}

   Here, :math:`n > p - 1` is the number of degrees of freedom, :math:`V` is a symmetric positive definite
   :math:`p`-by-:math:`p` scale matrix, whose Cholesky factor is specified by :data:`L`, and :data:`work` is
   :math:`p`-by-:math:`p` workspace. The :math:`p`-by-:math:`p` Wishart distributed matrix :math:`X` is stored
   in :data:`result` on output.

.. function:: int gsl_ran_wishart_pdf (const gsl_matrix * X, const gsl_matrix * L_X, const double n, const gsl_matrix * L, double * result, gsl_matrix * work)
              int gsl_ran_wishart_log_pdf (const gsl_matrix * X, const gsl_matrix * L_X, const double n, const gsl_matrix * L, double * result, gsl_matrix * work)

   These functions compute :math:`p(X)` or :math:`\log{p(X)}` for the :math:`p`-by-:math:`p` matrix
   :data:`X`, whose Cholesky factor is specified in :data:`L_X`. The degrees of freedom is
   given by :data:`n`, the Cholesky factor of the scale matrix :math:`V` is specified in :data:`L`,
   and :data:`work` is :math:`p`-by-:math:`p` workspace. The probably density value is returned
   in :data:`result`.

|newpage|

Shuffling and Sampling
======================

The following functions allow the shuffling and sampling of a set of
objects.  The algorithms rely on a random number generator as a source of
randomness and a poor quality generator can lead to correlations in the
output.  In particular it is important to avoid generators with a short
period.  For more information see Knuth, v2, 3rd ed, Section 3.4.2,
"Random Sampling and Shuffling".

.. function:: void gsl_ran_shuffle (const gsl_rng * r, void * base, size_t n, size_t size)

   This function randomly shuffles the order of :data:`n` objects, each of
   size :data:`size`, stored in the array :code:`base[0..n-1]`.  The
   output of the random number generator :data:`r` is used to produce the
   permutation.  The algorithm generates all possible :math:`n!`
   permutations with equal probability, assuming a perfect source of random
   numbers.

   The following code shows how to shuffle the numbers from 0 to 51::

      int a[52];

      for (i = 0; i < 52; i++)
        {
          a[i] = i;
        }

      gsl_ran_shuffle (r, a, 52, sizeof (int));

.. function:: int gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size)

   This function fills the array :code:`dest[k]` with :data:`k` objects taken
   randomly from the :data:`n` elements of the array
   :code:`src[0..n-1]`.  The objects are each of size :data:`size`.  The
   output of the random number generator :data:`r` is used to make the
   selection.  The algorithm ensures all possible samples are equally
   likely, assuming a perfect source of randomness.

   The objects are sampled **without** replacement, thus each object can
   only appear once in :data:`dest`.  It is required that :data:`k` be less
   than or equal to :data:`n`.  The objects in :data:`dest` will be in the
   same relative order as those in :data:`src`.  You will need to call
   :code:`gsl_ran_shuffle(r, dest, n, size)` if you want to randomize the
   order.

   The following code shows how to select a random sample of three unique
   numbers from the set 0 to 99::

      double a[3], b[100];

      for (i = 0; i < 100; i++)
        {
          b[i] = (double) i;
        }

      gsl_ran_choose (r, a, 3, b, 100, sizeof (double));

.. function:: void gsl_ran_sample (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size)

   This function is like :func:`gsl_ran_choose` but samples :data:`k` items
   from the original array of :data:`n` items :data:`src` with replacement, so
   the same object can appear more than once in the output sequence
   :data:`dest`.  There is no requirement that :data:`k` be less than :data:`n`
   in this case.

Examples
========

The following program demonstrates the use of a random number generator
to produce variates from a distribution.  It prints 10 samples from the
Poisson distribution with a mean of 3.

.. include:: examples/randpoisson.c
   :code:

If the library and header files are installed under :file:`/usr/local`
(the default location) then the program can be compiled with these
options::

  $ gcc -Wall demo.c -lgsl -lgslcblas -lm

Here is the output of the program,

.. include:: examples/randpoisson.txt
   :code:

The variates depend on the seed used by the generator.  The seed for the
default generator type :type:`gsl_rng_default` can be changed with the
:macro:`GSL_RNG_SEED` environment variable to produce a different stream
of variates::

  $ GSL_RNG_SEED=123 ./a.out

giving output

.. include:: examples/randpoisson2.txt
   :code:

The following program generates a random walk in two dimensions.

.. include:: examples/randwalk.c
   :code:

:numref:`fig_rand-walk` shows the output from the program.

.. _fig_rand-walk:

.. figure:: /images/random-walk.png
   :scale: 60%

   Four 10-step random walks from the origin.

The following program computes the upper and lower cumulative
distribution functions for the standard normal distribution at
:math:`x = 2`.

.. include:: examples/cdf.c
   :code:

Here is the output of the program,

.. include:: examples/cdf.txt
   :code:

References and Further Reading
==============================

For an encyclopaedic coverage of the subject readers are advised to
consult the book "Non-Uniform Random Variate Generation" by Luc
Devroye.  It covers every imaginable distribution and provides hundreds
of algorithms.

* Luc Devroye, "Non-Uniform Random Variate Generation",
  Springer-Verlag, ISBN 0-387-96305-7.  Available online at
  http://cg.scs.carleton.ca/~luc/rnbookindex.html.

The subject of random variate generation is also reviewed by Knuth, who
describes algorithms for all the major distributions.

* Donald E. Knuth, "The Art of Computer Programming: Seminumerical
  Algorithms" (Vol 2, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896842.

The Particle Data Group provides a short review of techniques for
generating distributions of random numbers in the "Monte Carlo"
section of its Annual Review of Particle Physics.

* Review of Particle Properties,
  R.M. Barnett et al., Physical Review D54, 1 (1996)
  http://pdg.lbl.gov/.

The Review of Particle Physics is available online in postscript and pdf
format.

An overview of methods used to compute cumulative distribution functions
can be found in *Statistical Computing* by W.J. Kennedy and
J.E. Gentle. Another general reference is *Elements of Statistical
Computing* by R.A. Thisted.

* William E. Kennedy and James E. Gentle, Statistical Computing (1980),
  Marcel Dekker, ISBN 0-8247-6898-1.

* Ronald A. Thisted, Elements of Statistical Computing (1988), 
  Chapman & Hall, ISBN 0-412-01371-1.

The cumulative distribution functions for the Gaussian distribution
are based on the following papers,

* Rational Chebyshev Approximations Using Linear Equations,
  W.J. Cody, W. Fraser, J.F. Hart. Numerische Mathematik 12, 242--251 (1968).

* Rational Chebyshev Approximations for the Error Function,
  W.J. Cody. Mathematics of Computation 23, n107, 631--637 (July 1969).
