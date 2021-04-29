GSL - GNU Scientific Library
============================

This is GSL, the GNU Scientific Library, a collection of numerical
routines for scientific computing.

GSL is free software, you can redistribute it and/or modify it under
the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

About this repository
=====================

This repository contains a copy of the latest stable version of GSL, 
with two main additions:

- a cmake build system
- optional inclusion of the AMPL GSL bindings, to access GSL functions from [AMPL](https://ampl.com/)

CMake build instructions
========================

Build with AMPL bindings
------------------------

To build the AMPL bindings: 
1) Make sure that the ASL submodule is initialized (in ```/ampl/thirdparty/asl```).
2) Create a build directory and move there:
   ```
   mkdir build
   cd build
   ```
3) Initialize the build files with your desired [generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html), for example:
   ```
   cmake .. -G"Unix Makefiles"
   ```

4) Build the resulting build files accordingly, for example:
   ```
   make .
   ```

Build without AMPL bindings 
---------------------------

If building GSL without the AMPL bindings:

1) Create a build directory and move there:
   ```
   mkdir build
   cd build
   ```
2) Initialize the build files with your desired [generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html), 
   defining the variable ```NO_AMPL_BINDINGS````
   ```
   cmake .. -G"Unix Makefiles" -DNO_AMPL_BINDINGS=1
   ```

3) Build the resulting build files accordingly, for example:
   ```
   make .
   ```


Other build parameters
----------------------

1) Defining GSL_DISABLE_TESTS skips the test generation, for example:
   ```
   cmake .. -DGSL_DISABLE_TESTS=1
   ```

2) There is experimental support to build only a subset of the GSL library;
   note that this is only supported when *not* building the AMPL bindings.
   To generate only some libraries of the GSL, define BUILDLIBS when calling cmake
   as a *comma separated* list of GSL directories. An example follows:

   ```
   cmake .. -DBUILDLIBS=ode-initval2,linalg -DNO_AMPL_BINDINGS=1
   ```

   All the dependencies are handled automatically; please note that the tests
   often have additional dependencies therefore, to fully benefit from this facility,
   disable them. For example:

    ```
   cmake .. -DBUILDLIBS=ode-initval2,linalg -DNO_AMPL_BINDINGS=1 -DGSL_DISABLE_TESTS=1
   ```


Availability
============

The current stable version of GSL is always available from ftp.gnu.org
in the directory /pub/gnu/gsl.

A list of mirror sites can be found at http://www.gnu.org/order/ftp.html

Installation
============

GSL follows the standard GNU installation procedure.  Please consult
the INSTALL file in this distribution for more detailed instructions.

For information about specific platforms and compilers see the
"Compilation Notes" section in the INSTALL file.

More information about GSL
==========================

The project homepage is http://www.gnu.org/software/gsl/

See the NEWS file for recent changes to the library.

The GSL Manual has been published and can be ordered from most
bookstores. The publication details are,

  GNU Scientific Library Reference Manual - Revised Second Edition, 
  M. Galassi et al, ISBN 0954161734 (620 pages, paperback).

The money raised from sales of the manual helps support the
development of GSL.

A Japanese translation of the reference manual is available from the
GSL website above (thanks to Daisuke TOMINAGA).

Reporting Bugs
==============

A list of known bugs can be found in the BUGS file.  Details of
compilation problems can be found in the INSTALL file.

If you find a bug which is not listed in these files please report it
to bug-gsl@gnu.org.

All bug reports should include:

       The version number of GSL, and where you obtained it.
       The hardware and operating system
       The compiler used, including version number and compilation options
       A description of the bug behaviour
       A short program which reproducibly exercises the bug

It is useful if you can check whether the same problem occurs when the
library is compiled without optimization.  Thank you.

Any errors or omissions in the manual can also be reported to the
same address.

Contributing to GSL
===================

If you are interested in participating in GSL development, please see
the webpage at http://www.gnu.org/software/gsl/

