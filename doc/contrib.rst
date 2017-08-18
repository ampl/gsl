*******************
Contributors to GSL
*******************

.. include:: include.rst

(See the :file:`AUTHORS` file in the distribution for up-to-date information.)

Mark Galassi

  Conceived GSL (with James Theiler) and wrote the design document.  Wrote
  the simulated annealing package and the relevant chapter in the manual.

James Theiler

  Conceived GSL (with Mark Galassi).  Wrote the random number generators
  and the relevant chapter in this manual.

Jim Davies

  Wrote the statistical routines and the relevant chapter in this
  manual.

Brian Gough

  FFTs, numerical integration, random number generators and distributions,
  root finding, minimization and fitting, polynomial solvers, complex
  numbers, physical constants, permutations, vector and matrix functions,
  histograms, statistics, ieee-utils, revised |cblas| Level 2 & 3,
  matrix decompositions, eigensystems, cumulative distribution functions,
  testing, documentation and releases.

Reid Priedhorsky

  Wrote and documented the initial version of the root finding routines
  while at Los Alamos National Laboratory, Mathematical Modeling and
  Analysis Group.  

.. email: reid@reidster.net

Gerard Jungman

  Special Functions, Series acceleration, ODEs, BLAS, Linear Algebra,
  Eigensystems, Hankel Transforms.

Patrick Alken

  Implementation of nonsymmetric and generalized eigensystems, B-splines,
  linear and nonlinear least squares, matrix decompositions,
  associated Legendre functions, running statistics, sparse matrices,
  and sparse linear algebra.

Mike Booth

  Wrote the Monte Carlo library.

Jorma Olavi Tähtinen

  Wrote the initial complex arithmetic functions.

Thomas Walter 

  Wrote the initial heapsort routines and Cholesky decomposition.

Fabrice Rossi

  Multidimensional minimization.

Carlo Perassi

  Implementation of the random number generators in Knuth's
  *Seminumerical Algorithms*, 3rd Ed.

Szymon Jaroszewicz 

  Wrote the routines for generating combinations.

.. <sj@cs.umb.edu>

Nicolas Darnis

  Wrote the cyclic functions and the initial functions for canonical 
  permutations.

Jason H. Stover

  Wrote the major cumulative distribution functions.

.. (jason@sakla.net) 

Ivo Alxneit

  Wrote the routines for wavelet transforms.

Tuomo Keskitalo

  Improved the implementation of the ODE solvers and wrote the
  ode-initval2 routines.

Lowell Johnson

  Implementation of the Mathieu functions.

Rhys Ulerich 

  Wrote the multiset routines.

Pavel Holoborodko 

  Wrote the fixed order Gauss-Legendre quadrature routines.

Pedro Gonnet

  Wrote the |cquad| integration routines.

Thanks to Nigel Lowry for help in proofreading the manual.

The non-symmetric eigensystems routines contain code based on the
LAPACK linear algebra library.  LAPACK is distributed under the
following license::

  Copyright (c) 1992-2006 The University of Tennessee.  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer. 
  
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer listed
    in this license in the documentation and/or other materials
    provided with the distribution.
  
  * Neither the name of the copyright holders nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
