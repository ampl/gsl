.. index:: Mathieu functions

The routines described in this section compute the angular and radial
Mathieu functions, and their characteristic values.  Mathieu
functions are the solutions of the following two differential
equations:

.. only:: not texinfo

   .. math::

      {{d^2 y}\over{d v^2}}& + (a - 2q\cos 2v)y  = 0 \\
      {{d^2 f}\over{d u^2}}& - (a - 2q\cosh 2u)f  = 0

.. only:: texinfo

   ::

      d^2y/dv^2 + (a - 2q\cos 2v)y = 0
      d^2f/du^2 - (a - 2q\cosh 2u)f = 0

The angular Mathieu functions :math:`ce_r(x,q)`, :math:`se_r(x,q)` are
the even and odd periodic solutions of the first equation, which is known as Mathieu's equation. These exist
only for the discrete sequence of  characteristic values :math:`a = a_r(q)`
(even-periodic) and :math:`a = b_r(q)` (odd-periodic).

The radial Mathieu functions :math:`Mc^{(j)}_{r}(z,q)` and
:math:`Ms^{(j)}_{r}(z,q)`
are the solutions of the second equation,
which is referred to as Mathieu's modified equation.  The
radial Mathieu functions of the first, second, third and fourth kind
are denoted by the parameter :math:`j`, which takes the value 1, 2, 3
or 4.

.. The angular Mathieu functions can be divided into four types as
.. @tex
.. \beforedisplay
.. $$
.. \eqalign{
.. x & = \sum_{m=0}^\infty A_{2m+p} \cos(2m+p)\phi, \quad p = 0, 1, \cr
.. x & = \sum_{m=0}^\infty B_{2m+p} \sin(2m+p)\phi, \quad p = 0, 1.
.. }
.. $$
.. \afterdisplay
.. @end tex
.. @ifinfo

.. @example
.. x = \sum_(m=0)^\infty A_(2m+p) \cos(2m+p)\phi,   p = 0, 1,
.. x = \sum_(m=0)^\infty B_(2m+p) \sin(2m+p)\phi,   p = 0, 1.
.. @end example

.. @end ifinfo
.. @noindent
.. The nomenclature used for the angular Mathieu functions is :math:`ce_n`
.. for the first solution and :math:`se_n` for the second.

.. Similar solutions exist for the radial Mathieu functions by replacing
.. the trigonometric functions with their corresponding hyperbolic
.. functions as shown below.
.. @tex
.. \beforedisplay
.. $$
.. \eqalign{
.. x & = \sum_{m=0}^\infty A_{2m+p} \cosh(2m+p)u, \quad p = 0, 1, \cr
.. x & = \sum_{m=0}^\infty B_{2m+p} \sinh(2m+p)u, \quad p = 0, 1.
.. }
.. $$
.. \afterdisplay
.. @end tex
.. @ifinfo

.. @example
.. x = \sum_(m=0)^\infty A_(2m+p) \cosh(2m+p)u,   p = 0, 1,
.. x = \sum_(m=0)^\infty B_(2m+p) \sinh(2m+p)u,   p = 0, 1.
.. @end example

.. @end ifinfo
.. @noindent
.. The nomenclature used for the radial Mathieu functions is :math:`Mc_n`
.. for the first solution and :math:`Ms_n` for the second.  The hyperbolic
.. series do not always converge at an acceptable rate.  Therefore most
.. texts on the subject suggest using the following equivalent equations
.. that are expanded in series of Bessel and Hankel functions.
.. @tex
.. \beforedisplay
.. $$
.. \eqalign{
.. Mc_{2n}^{(j)}(x,q) & = \sum_{m=0}^\infty (-1)^{r+k}
..       A_{2m}^{2n}(q)\left[J_m(u_1)Z_m^{(j)}(u_2) +
..                           J_m(u_1)Z_m^{(j)}(u_2)\right]/A_2^{2n} \cr
.. Mc_{2n+1}^{(j)}(x,q) & = \sum_{m=0}^\infty (-1)^{r+k}
..       A_{2m+1}^{2n+1}(q)\left[J_m(u_1)Z_{m+1}^{(j)}(u_2) +
..                               J_{m+1}(u_1)Z_m^{(j)}(u_2)\right]/A_1^{2n+1} \cr
.. Ms_{2n}^{(j)}(x,q) & = \sum_{m=1}^\infty (-1)^{r+k}
..       B_{2m}^{2n}(q)\left[J_{m-1}(u_1)Z_{m+1}^{(j)}(u_2) +
..                           J_{m+1}(u_1)Z_{m-1}^{(j)}(u_2)\right]/B_2^{2n} \cr
.. Ms_{2n+1}^{(j)}(x,q) & = \sum_{m=0}^\infty (-1)^{r+k}
..       B_{2m+1}^{2n+1}(q)\left[J_m(u_1)Z_{m+1}^{(j)}(u_2) +
..                               J_{m+1}(u_1)Z_m^{(j)}(u_2)\right]/B_1^{2n+1}
.. }
.. $$
.. \afterdisplay
.. @end tex
.. @ifinfo

.. @example
.. Mc_(2n)^(j)(x,q) = \sum_(m=0)^\infty (-1)^(r+k) A_(2m)^(2n)(q)
..     [J_m(u_1)Z_m^(j)(u_2) + J_m(u_1)Z_m^(j)(u_2)]/A_2^(2n)
.. Mc_(2n+1)^(j)(x,q) = \sum_(m=0)^\infty (-1)^(r+k) A_(2m+1)^(2n+1)(q)
..     [J_m(u_1)Z_(m+1)^(j)(u_2) + J_(m+1)(u_1)Z_m^(j)(u_2)]/A_1^(2n+1)
.. Ms_(2n)^(j)(x,q) = \sum_(m=1)^\infty (-1)^(r+k) B_(2m)^(2n)(q)
..     [J_(m-1)(u_1)Z_(m+1)^(j)(u_2) + J_(m+1)(u_1)Z_(m-1)^(j)(u_2)]/B_2^(2n)
.. Ms_(2n+1)^(j)(x,q) = \sum_(m=0)^\infty (-1)^(r+k) B_(2m+1)^(2n+1)(q)
..     [J_m(u_1)Z_(m+1)^(j)(u_2) + J_(m+1)(u_1)Z_m^(j)(u_2)]/B_1^(2n+1)
.. @end example

.. @end ifinfo
.. @noindent
.. where @c{$u_1 = \sqrt{q} \exp(-x)$} 
.. @math{u_1 = \sqrt@{q@} \exp(-x)} and @c{$u_2 = \sqrt@{q@} \exp(x)$}
.. @math{u_2 = \sqrt@{q@} \exp(x)} and
.. @tex
.. \beforedisplay
.. $$
.. \eqalign{
.. Z_m^{(1)}(u) & = J_m(u) \cr
.. Z_m^{(2)}(u) & = Y_m(u) \cr
.. Z_m^{(3)}(u) & = H_m^{(1)}(u) \cr
.. Z_m^{(4)}(u) & = H_m^{(2)}(u)
.. }
.. $$
.. \afterdisplay
.. @end tex
.. @ifinfo

.. @example
.. Z_m^(1)(u) = J_m(u)
.. Z_m^(2)(u) = Y_m(u)
.. Z_m^(3)(u) = H_m^(1)(u)
.. Z_m^(4)(u) = H_m^(2)(u)
.. @end example

.. @end ifinfo
.. @noindent
.. where @math{J_m(u)}, @math{Y_m(u)}, @math{H_m^{(1)}(u)}, and
.. :math:`H_m^{(2)}(u)` are the regular and irregular Bessel functions and
.. the Hankel functions, respectively.

For more information on the Mathieu functions, see Abramowitz and
Stegun, Chapter 20.  These functions are defined in the header file
:file:`gsl_sf_mathieu.h`.

Mathieu Function Workspace
--------------------------

The Mathieu functions can be computed for a single order or
for multiple orders, using array-based routines.  The array-based
routines require a preallocated workspace.

.. type:: gsl_sf_mathieu_workspace

   Workspace required for array-based routines

.. function:: gsl_sf_mathieu_workspace * gsl_sf_mathieu_alloc (size_t n, double qmax)

   This function returns a workspace for the array versions of the
   Mathieu routines.  The arguments n and :data:`qmax` specify the
   maximum order and :math:`q`-value of Mathieu functions which can be
   computed with this workspace.  

.. This is required in order to properly
.. terminate the infinite eigenvalue matrix for high precision solutions.
.. The characteristic values for all orders :math:`0 \to n` are stored in
.. the work structure array element @kbd{work->char_value}.

.. function:: void gsl_sf_mathieu_free (gsl_sf_mathieu_workspace * work)

   This function frees the workspace :data:`work`.

Mathieu Function Characteristic Values
--------------------------------------
.. index:: Mathieu Function Characteristic Values

.. function:: int gsl_sf_mathieu_a (int n, double q)
              int gsl_sf_mathieu_a_e (int n, double q, gsl_sf_result * result)
              int gsl_sf_mathieu_b (int n, double q)
              int gsl_sf_mathieu_b_e (int n, double q, gsl_sf_result * result)

   These routines compute the characteristic values :math:`a_n(q)`,
   :math:`b_n(q)` of the Mathieu functions :math:`ce_n(q,x)` and
   :math:`se_n(q,x)`, respectively.

.. function:: int gsl_sf_mathieu_a_array (int order_min, int order_max, double q, gsl_sf_mathieu_workspace * work, double result_array[])
              int gsl_sf_mathieu_b_array (int order_min, int order_max, double q, gsl_sf_mathieu_workspace * work, double result_array[])

   These routines compute a series of Mathieu characteristic values
   :math:`a_n(q)`, :math:`b_n(q)` for :math:`n` from :data:`order_min` to
   :data:`order_max` inclusive, storing the results in the array :data:`result_array`.

Angular Mathieu Functions
-------------------------
.. index::
   single: Angular Mathieu Functions
   single: ce(q,x), Mathieu function
   single: se(q,x), Mathieu function

.. function:: int gsl_sf_mathieu_ce (int n, double q, double x)
              int gsl_sf_mathieu_ce_e (int n, double q, double x, gsl_sf_result * result)
              int gsl_sf_mathieu_se (int n, double q, double x)
              int gsl_sf_mathieu_se_e (int n, double q, double x, gsl_sf_result * result)

   These routines compute the angular Mathieu functions :math:`ce_n(q,x)`
   and :math:`se_n(q,x)`, respectively.

.. function:: int gsl_sf_mathieu_ce_array (int nmin, int nmax, double q, double x, gsl_sf_mathieu_workspace * work, double result_array[])
              int gsl_sf_mathieu_se_array (int nmin, int nmax, double q, double x, gsl_sf_mathieu_workspace * work, double result_array[])

   These routines compute a series of the angular Mathieu functions
   :math:`ce_n(q,x)` and :math:`se_n(q,x)` of order :math:`n` from
   :data:`nmin` to :data:`nmax` inclusive, storing the results in the array
   :data:`result_array`.

Radial Mathieu Functions
------------------------
.. index:: Radial Mathieu Functions

.. function:: int gsl_sf_mathieu_Mc (int j, int n, double q, double x)
              int gsl_sf_mathieu_Mc_e (int j, int n, double q, double x, gsl_sf_result * result)
              int gsl_sf_mathieu_Ms (int j, int n, double q, double x)
              int gsl_sf_mathieu_Ms_e (int j, int n, double q, double x, gsl_sf_result * result)

   These routines compute the radial :data:`j`-th kind Mathieu functions
   :math:`Mc_n^{(j)}(q,x)` and :math:`Ms_n^{(j)}(q,x)` of order :data:`n`.

   The allowed values of :data:`j` are 1 and 2.
   The functions for :math:`j = 3,4` can be computed as 
   :math:`M_n^{(3)} = M_n^{(1)} + iM_n^{(2)}` and
   :math:`M_n^{(4)} = M_n^{(1)} - iM_n^{(2)}`,
   where 
   :math:`M_n^{(j)} = Mc_n^{(j)}` or
   :math:`Ms_n^{(j)}`.

.. function:: int gsl_sf_mathieu_Mc_array (int j, int nmin, int nmax, double q, double x, gsl_sf_mathieu_workspace * work, double result_array[])
              int gsl_sf_mathieu_Ms_array (int j, int nmin, int nmax, double q, double x, gsl_sf_mathieu_workspace * work, double result_array[])

   These routines compute a series of the radial Mathieu functions of
   kind :data:`j`, with order from :data:`nmin` to :data:`nmax` inclusive, storing the
   results in the array :data:`result_array`.
