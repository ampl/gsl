.. index::
   license of GSL
   GNU General Public License

************
Introduction
************

The GNU Scientific Library (GSL) is a collection of routines for
numerical computing.  The routines have been written from scratch in C,
and present a modern Applications Programming Interface
(API) for C programmers, allowing wrappers to be written for very
high level languages.  The source code is distributed under the GNU
General Public License.

Routines available in GSL
=========================

The library covers a wide range of topics in numerical computing.
Routines are available for the following areas,

===========================  ===========================  ===========================
Complex Numbers              Roots of Polynomials         Special Functions
Vectors and Matrices         Permutations                 Combinations
Sorting                      BLAS Support                 Linear Algebra
CBLAS Library                Fast Fourier Transforms      Eigensystems
Random Numbers               Quadrature                   Random Distributions
Quasi-Random Sequences       Histograms                   Statistics
Monte Carlo Integration      N-Tuples                     Differential Equations
Simulated Annealing          Numerical Differentiation    Interpolation
Series Acceleration          Chebyshev Approximations     Root-Finding
Discrete Hankel Transforms   Least-Squares Fitting        Minimization
IEEE Floating-Point          Physical Constants           Basis Splines
Wavelets                     Sparse BLAS Support          Sparse Linear Algebra
===========================  ===========================  ===========================

The use of these routines is described in this manual.  Each chapter
provides detailed definitions of the functions, followed by example
programs and references to the articles on which the algorithms are
based.

Where possible the routines have been based on reliable public-domain
packages such as FFTPACK and QUADPACK, which the developers of GSL
have reimplemented in C with modern coding conventions.

.. index::
   single: free software, explanation of

GSL is Free Software
====================

The subroutines in the GNU Scientific Library are "free software";
this means that everyone is free to use them, and to redistribute them
in other free programs.  The library is not in the public domain; it is
copyrighted and there are conditions on its distribution.  These
conditions are designed to permit everything that a good cooperating
citizen would want to do.  What is not allowed is to try to prevent
others from further sharing any version of the software that they might
get from you.

Specifically, we want to make sure that you have the right to share
copies of programs that you are given which use the GNU Scientific
Library, that you receive their source code or else can get it if you
want it, that you can change these programs or use pieces of them in new
free programs, and that you know you can do these things.

To make sure that everyone has such rights, we have to forbid you to
deprive anyone else of these rights.  For example, if you distribute
copies of any code which uses the GNU Scientific Library, you must give
the recipients all the rights that you have received.  You must make
sure that they, too, receive or can get the source code, both to the
library and the code which uses it.  And you must tell them their
rights.  This means that the library should not be redistributed in
proprietary programs.

Also, for our own protection, we must make certain that everyone finds
out that there is no warranty for the GNU Scientific Library.  If these
programs are modified by someone else and passed on, we want their
recipients to know that what they have is not what we distributed, so
that any problems introduced by others will not reflect on our
reputation.

The precise conditions for the distribution of software related to the
GNU Scientific Library are found in the
`GNU General Public License <https://www.gnu.org/software/gsl/manual/html_node/GNU-General-Public-License.html#GNU-General-Public-License>`_.
Further information about this
license is available from the GNU Project webpage `Frequently Asked
Questions about the GNU GPL <http://www.gnu.org/copyleft/gpl-faq.html>`_.

The Free Software Foundation also operates a license consulting
service for commercial users (contact details available from
http://www.fsf.org.

.. index::
   obtaining GSL
   downloading GSL
   mailing list for GSL announcements
   info-gsl mailing list

Obtaining GSL
=============

The source code for the library can be obtained in different ways, by
copying it from a friend, purchasing it on CDROM or downloading it
from the internet. A list of public ftp servers which carry the source
code can be found on the GNU website, http://www.gnu.org/software/gsl/.

The preferred platform for the library is a GNU system, which allows it
to take advantage of additional features in the GNU C compiler and GNU C
library.  However, the library is fully portable and should compile on
most systems with a C compiler. 

Announcements of new releases, updates and other relevant events are
made on the info-gsl@gnu.org mailing list.  To subscribe to this
low-volume list, send an email of the following form::

    To: info-gsl-request@gnu.org 
    Subject: subscribe

You will receive a response asking you to reply in order to confirm
your subscription.

.. index::
   warranty (none)

No Warranty
===========

The software described in this manual has no warranty, it is provided
"as is".  It is your responsibility to validate the behavior of the
routines and their accuracy using the source code provided, or to
purchase support and warranties from commercial redistributors.  Consult the
`GNU General Public License <https://www.gnu.org/software/gsl/manual/html_node/GNU-General-Public-License.html#GNU-General-Public-License>`_
for further details.

.. index::
   reporting bugs in GSL
   bugs, how to report
   bug-gsl mailing list
   mailing list, bug-gsl

Reporting Bugs
==============

A list of known bugs can be found in the :file:`BUGS` file included in
the GSL distribution or online in the GSL bug tracker. [#f1]_
Details of compilation problems can be found in the :file:`INSTALL` file.

If you find a bug which is not listed in these files, please report it to
bug-gsl@gnu.org.

All bug reports should include:

- The version number of GSL
- The hardware and operating system
- The compiler used, including version number and compilation options
- A description of the bug behavior
- A short program which exercises the bug

It is useful if you can check whether the same problem occurs when the
library is compiled without optimization.  Thank you.

Any errors or omissions in this manual can also be reported to the
same address.

.. index::
   mailing list archives
   single: website, developer information
   contacting the GSL developers

Further Information
===================

Additional information, including online copies of this manual, links to
related projects, and mailing list archives are available from the
website mentioned above.  

Any questions about the use and installation of the library can be asked
on the mailing list help-gsl@gnu.org.  To subscribe to this
list, send an email of the following form::

    To: help-gsl-request@gnu.org
    Subject: subscribe

This mailing list can be used to ask questions not covered by this
manual, and to contact the developers of the library.

If you would like to refer to the GNU Scientific Library in a journal
article, the recommended way is to cite this reference manual,
e.g.::

    M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078.

If you want to give a url, use "http://www.gnu.org/software/gsl/".

.. index::
   single: conventions, used in manual
   single: examples, conventions used in
   single: shell prompt
   single: $, shell prompt

Conventions used in this manual
===============================

.. index::
   single: dollar sign $, shell prompt

This manual contains many examples which can be typed at the keyboard.
A command entered at the terminal is shown like this::

    $ command

The first character on the line is the terminal prompt, and should not
be typed.  The dollar sign $ is used as the standard prompt in
this manual, although some systems may use a different character.

The examples assume the use of the GNU operating system.  There may be
minor differences in the output on other systems.  The commands for
setting environment variables use the Bourne shell syntax of the
standard GNU shell (:code:`bash`).

.. rubric:: Footnotes

.. [#f1] http://savannah.gnu.org/bugs/?group=gsl
