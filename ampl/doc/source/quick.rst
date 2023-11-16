Quick Start
===========

This library extends AMPL with over 300 functions chosen from the `GNU
Scientific Library <http://www.gnu.org/software/gsl/>`__, including
scientific, statistical, and utility functions of many kinds. 
Once installed and declared, these functions may be applied to any arguments
anywhere in an AMPL model or script.

When applied to variables in objectives or constraints, these functions
are included in AMPL’s nonlinear expression representations, and are
accessible to solvers. Thus solvers that rely on function evaluations
can optimize over expressions that use this library. This feature is
most useful for applying local nonlinear solvers — such as CONOPT,
Knitro, LOQO, MINOS, SNOPT — to continuous functions of one variable.
Many of the library’s continuous functions also provide the 1st and 2nd
derivative evaluations that some solvers require.

For over 30 probability distributions, the library provides random
variates, probability densities, and cumulative distributions and their
inverses. Functions to compute means, variances, correlations, and other
statistical functions are also provided. A sort function can be used to
rearrange the values in an indexed AMPL parameter or variable.

Downloading and installing the library
--------------------------------------

The function library is packaged in a single file named ``amplgsl.dll``.
This file is available for Windows, Linux, and macOS 
from `portal.ampl.com <https://portal.ampl.com/>`__ in 
`My Downloads <https://portal.ampl.com/user/ampl/download/list>`__.

To install, use one of the above links to download the appropriate
``tgz`` or ``zip`` file, and copy that file to the same folder or
directory as your AMPL executable (``ampl`` or ``ampl.exe``). Then
extract ``amplgsl.dll`` by double-clicking the downloaded file’s icon on
Windows and macOS systems, or by invoking

.. code:: bash

   gzip -dc amplgsl.linux-intel64.tgz | tar xf -

in a Linux command window.

Using the functions in AMPL
---------------------------

Once in each session that uses the extended function library, the
library file must be loaded with the command

.. code:: ampl

   load amplgsl.dll;

Subsequently, before a particular function is used, it must be declared
to the current model with a command of the form

.. code:: ampl

   function gsl_ _funcname_ ;

where ``_funcname_`` is replaced by any of the names in the 
:doc:`library listing <reference>`.
As an example, the following model, based on the Hock-Schittkowski test
set, uses a unit Gaussian cumulative distribution function from the
extended library:

.. code:: ampl

   function gsl_cdf_ugaussian_P;  

   param l {1..4};
   param u {1..4};
   var x {j in 1..4} >= l[j], <= u[j] := 1;  

   param a = 0.1;
   param b = 1000;
   param d = 1;
   param n = 4;  

   minimize obj:
     ( a*n - (b*(exp(x[1])-1) - x[3])*x[4]/(exp(x[1]) - 1 + x[4]) )/x[1] ;  

   subject to constr1:
       x[3] = 2*gsl_cdf_ugaussian_P(-x[2]);  

   subject to constr2:
       x[4] = gsl_cdf_ugaussian_P(-x[2] + d*sqrt(n)) +
              gsl_cdf_ugaussian_P(-x[2] - d*sqrt(n)); 

With the appropriate data file,

.. code:: ampl

   param:    l       u :=
     1    0.0001   100
     2    0        100
     3    0          2
     4    0          2 ;

the problem can be solved like this:

.. code:: ampl

   ampl: model hs069.mod;
   ampl: data hs069.dat;
   ampl: load amplgsl.dll;  

   ampl: option solver knitro;
   ampl: solve;
   Knitro 12.4.0: Locally optimal or satisfactory solution.
   objective -956.7128869; feasibility error 5.49e-09
   10 iterations; 11 function evaluations  

   ampl: display x;
   1  0.0293714
   2  1.19025
   3  0.233947
   4  0.791668

As a convenience, ``function`` statements for all members of the
Extended Function Library are provided as part of the library
distribution, in a file named ``gsl.ampl``. Thus the single command
``include gsl.ampl;`` makes all of the library functions available for
use in an AMPL session.

Documentation and licensing
---------------------------

Complete documentation for this library is provided in 
:doc:`AMPL Bindings for the GNU Scientific Library <reference>`. 
It includes a listing of all functions and their definitions, 
further instructions for use, and complete license details.

This free-standing library is provided under the `GNU General Public
License <http://www.gnu.org/licenses/gpl.html>`__ (GPL). In general
terms, the GNU GPL does not restrict the use of this library for
internal purposes within any organization (whether a company,
university, or government agency). However if this library is included
within a “release” of software to the public, then you may be required
to offer that release under the GPL. See the extensive list of
`Frequently Asked Questions about the GNU
Licenses <http://www.gnu.org/licenses/gpl-faq.html>`__ for more
information about these provisions.

For uses that may require functions covered by more permissive
open-source licenses, or functions from commercial libraries, please
contact us at support@ampl.com for assistance with alternatives.
