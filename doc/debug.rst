****************************
Debugging Numerical Programs
****************************

This chapter describes some tips and tricks for debugging numerical
programs which use GSL.

.. index::
   single: gdb
   single: debugging numerical programs
   single: breakpoints

Using gdb
=========

Any errors reported by the library are passed to the function
:func:`gsl_error`.  By running your programs under gdb and setting a
breakpoint in this function you can automatically catch any library
errors.  You can add a breakpoint for every session by putting::

  break gsl_error

into your :file:`.gdbinit` file in the directory where your program is
started.  

If the breakpoint catches an error then you can use a backtrace
(:code:`bt`) to see the call-tree, and the arguments which possibly
caused the error.  By moving up into the calling function you can
investigate the values of variables at that point.  Here is an example
from the program :code:`fft/test_trap`, which contains the following
line::

  status = gsl_fft_complex_wavetable_alloc (0, &complex_wavetable);

The function :func:`gsl_fft_complex_wavetable_alloc` takes the length of
an FFT as its first argument.  When this line is executed an error will
be generated because the length of an FFT is not allowed to be zero.

To debug this problem we start :code:`gdb`, using the file
:file:`.gdbinit` to define a breakpoint in :func:`gsl_error`::

  $ gdb test_trap

  GDB is free software and you are welcome to distribute copies
  of it under certain conditions; type "show copying" to see
  the conditions.  There is absolutely no warranty for GDB;
  type "show warranty" for details.  GDB 4.16 (i586-debian-linux), 
  Copyright 1996 Free Software Foundation, Inc.

  Breakpoint 1 at 0x8050b1e: file error.c, line 14.

When we run the program this breakpoint catches the error and shows the
reason for it::

  (gdb) run
  Starting program: test_trap 

  Breakpoint 1, gsl_error (reason=0x8052b0d 
      "length n must be positive integer", 
      file=0x8052b04 "c_init.c", line=108, gsl_errno=1) 
      at error.c:14
  14        if (gsl_error_handler) 

The first argument of :func:`gsl_error` is always a string describing the
error.  Now we can look at the backtrace to see what caused the problem::

  (gdb) bt
  #0  gsl_error (reason=0x8052b0d 
      "length n must be positive integer", 
      file=0x8052b04 "c_init.c", line=108, gsl_errno=1)
      at error.c:14
  #1  0x8049376 in gsl_fft_complex_wavetable_alloc (n=0,
      wavetable=0xbffff778) at c_init.c:108
  #2  0x8048a00 in main (argc=1, argv=0xbffff9bc) 
      at test_trap.c:94
  #3  0x80488be in ___crt_dummy__ ()

We can see that the error was generated in the function
:func:`gsl_fft_complex_wavetable_alloc` when it was called with an
argument of :code:`n = 0`.  The original call came from line 94 in the
file :file:`test_trap.c`.

By moving up to the level of the original call we can find the line that
caused the error::

  (gdb) up
  #1  0x8049376 in gsl_fft_complex_wavetable_alloc (n=0,
      wavetable=0xbffff778) at c_init.c:108
  108   GSL_ERROR ("length n must be positive integer", GSL_EDOM);
  (gdb) up
  #2  0x8048a00 in main (argc=1, argv=0xbffff9bc) 
    at test_trap.c:94
  94    status = gsl_fft_complex_wavetable_alloc (0,
          &complex_wavetable);

Thus we have found the line that caused the problem.  From this point we
could also print out the values of other variables such as
:code:`complex_wavetable`.

.. index:: floating point registers

Examining floating point registers
==================================

The contents of floating point registers can be examined using the
command :code:`info float` (on supported platforms)::

  (gdb) info float
       st0: 0xc4018b895aa17a945000  Valid Normal -7.838871e+308
       st1: 0x3ff9ea3f50e4d7275000  Valid Normal 0.0285946
       st2: 0x3fe790c64ce27dad4800  Valid Normal 6.7415931e-08
       st3: 0x3ffaa3ef0df6607d7800  Spec  Normal 0.0400229
       st4: 0x3c028000000000000000  Valid Normal 4.4501477e-308
       st5: 0x3ffef5412c22219d9000  Zero  Normal 0.9580257
       st6: 0x3fff8000000000000000  Valid Normal 1
       st7: 0xc4028b65a1f6d243c800  Valid Normal -1.566206e+309
     fctrl: 0x0272 53 bit; NEAR; mask DENOR UNDER LOS;
     fstat: 0xb9ba flags 0001; top 7; excep DENOR OVERF UNDER LOS
      ftag: 0x3fff
       fip: 0x08048b5c
       fcs: 0x051a0023
    fopoff: 0x08086820
    fopsel: 0x002b

Individual registers can be examined using the variables :code:`$reg`,
where :code:`reg` is the register name::

  (gdb) p $st1 
  $1 = 0.02859464454261210347719

.. index::
   single: exceptions, floating point
   single: floating point exceptions

Handling floating point exceptions
==================================

It is possible to stop the program whenever a :code:`SIGFPE` floating
point exception occurs.  This can be useful for finding the cause of an
unexpected infinity or :code:`NaN`.  The current handler settings can be
shown with the command :code:`info signal SIGFPE`::

  (gdb) info signal SIGFPE
  Signal  Stop  Print  Pass to program Description
  SIGFPE  Yes   Yes    Yes             Arithmetic exception

Unless the program uses a signal handler the default setting should be
changed so that SIGFPE is not passed to the program, as this would cause
it to exit.  The command :code:`handle SIGFPE stop nopass` prevents this::

  (gdb) handle SIGFPE stop nopass
  Signal  Stop  Print  Pass to program Description
  SIGFPE  Yes   Yes    No              Arithmetic exception

Depending on the platform it may be necessary to instruct the kernel to
generate signals for floating point exceptions.  For programs using GSL
this can be achieved using the :macro:`GSL_IEEE_MODE` environment variable
in conjunction with the function :func:`gsl_ieee_env_setup` as described
in :ref:`chap_ieee`::

  (gdb) set env GSL_IEEE_MODE=double-precision

.. index::
   single: warning options
   single: gcc warning options

GCC warning options for numerical programs
==========================================

Writing reliable numerical programs in C requires great care.  The
following GCC warning options are recommended when compiling numerical
programs::

  gcc -ansi -pedantic -Werror -Wall -W 
    -Wmissing-prototypes -Wstrict-prototypes 
    -Wconversion -Wshadow -Wpointer-arith 
    -Wcast-qual -Wcast-align 
    -Wwrite-strings -Wnested-externs 
    -fshort-enums -fno-common -Dinline= -g -O2

.. Uninitialized variables, conversions to and from integers or from
.. signed to unsigned integers can all cause hard-to-find problems.  For
.. many non-numerical programs compiling with :code:`gcc`'s warning option
.. :code:`-Wall` provides a good check against common errors.  However, for
.. numerical programs :code:`-Wall` is not enough. 

.. If you are unconvinced take a look at this program which contains an
.. error that can occur in numerical code,

.. @example
.. #include <math.h>
.. #include <stdio.h>

.. double f (int x);

.. int main ()
.. @{
..   double a = 1.5;
..   double y = f(a);
..   printf("a = %g, sqrt(a) = %g\n", a, y);  
..   return 0;
.. @}

.. double f(x) @{
..   return sqrt(x);
.. @}
.. @end example

.. @noindent
.. This code compiles cleanly with :code:`-Wall` but produces some strange
.. output,

.. @example
.. bash$ gcc -Wall tmp.c -lm
.. bash$ ./a.out 
.. a = 1.5, sqrt(a) = 1
.. @end example

.. @noindent
.. Note that adding :code:`-ansi` does not help here, since the program does
.. not contain any invalid constructs.  What is happening is that the
.. prototype for the function :code:`f(int x)` is not consistent with the
.. function call :code:`f(y)`, where :code:`y` is a floating point
.. number.  This results in the argument being silently converted to an
.. integer.  This is valid C, but in a numerical program it also likely to
.. be a programming error so we would like to be warned about it. (If we
.. genuinely wanted to convert :code:`y` to an integer then we could use an
.. explicit cast, :code:`(int)y`).  

.. Fortunately GCC provides many additional warnings which can alert you to
.. problems such as this.  You just have to remember to use them.  Here is a
.. set of recommended warning options for numerical programs.

For details of each option consult the manual *Using and Porting
GCC*.  The following table gives a brief explanation of what types of
errors these options catch.

:code:`-ansi -pedantic`

  Use ANSI C, and reject any non-ANSI extensions.  These flags help in
  writing portable programs that will compile on other systems.

:code:`-Werror`

  Consider warnings to be errors, so that compilation stops.  This prevents
  warnings from scrolling off the top of the screen and being lost.  You
  won't be able to compile the program until it is completely
  warning-free.

:code:`-Wall`

  This turns on a set of warnings for common programming problems.  You
  need :code:`-Wall`, but it is not enough on its own.

:code:`-O2`

  Turn on optimization.  The warnings for uninitialized variables in
  :code:`-Wall` rely on the optimizer to analyze the code.  If there is no
  optimization then these warnings aren't generated.

:code:`-W`

  This turns on some extra warnings not included in :code:`-Wall`, such as
  missing return values and comparisons between signed and unsigned
  integers.

:code:`-Wmissing-prototypes -Wstrict-prototypes`

  Warn if there are any missing or inconsistent prototypes.  Without
  prototypes it is harder to detect problems with incorrect arguments.

:code:`-Wconversion`

  The main use of this option is to warn about conversions from signed to
  unsigned integers.  For example, :code:`unsigned int x = -1`.  If you need
  to perform such a conversion you can use an explicit cast.

:code:`-Wshadow`

  This warns whenever a local variable shadows another local variable.  If
  two variables have the same name then it is a potential source of
  confusion.

:code:`-Wpointer-arith -Wcast-qual -Wcast-align`

  These options warn if you try to do pointer arithmetic for types which
  don't have a size, such as :code:`void`, if you remove a :code:`const`
  cast from a pointer, or if you cast a pointer to a type which has a
  different size, causing an invalid alignment.

:code:`-Wwrite-strings`

  This option gives string constants a :code:`const` qualifier so that it
  will be a compile-time error to attempt to overwrite them.

:code:`-fshort-enums`

  This option makes the type of :code:`enum` as short as possible.  Normally
  this makes an :code:`enum` different from an :code:`int`.  Consequently any
  attempts to assign a pointer-to-int to a pointer-to-enum will generate a
  cast-alignment warning.

:code:`-fno-common`

  This option prevents global variables being simultaneously defined in
  different object files (you get an error at link time).  Such a variable
  should be defined in one file and referred to in other files with an
  :code:`extern` declaration.

:code:`-Wnested-externs`

  This warns if an :code:`extern` declaration is encountered within a
  function.

:code:`-Dinline=`

  The :code:`inline` keyword is not part of ANSI C. Thus if you want to use
  :code:`-ansi` with a program which uses inline functions you can use this
  preprocessor definition to remove the :code:`inline` keywords.

:code:`-g`

  It always makes sense to put debugging symbols in the executable so that
  you can debug it using :code:`gdb`.  The only effect of debugging symbols
  is to increase the size of the file, and you can use the :code:`strip`
  command to remove them later if necessary.

.. For comparison, this is what happens when the test program above is
.. compiled with these options.

.. @example
.. bash$ gcc -ansi -pedantic -Werror -W -Wall -Wtraditional 
.. -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align 
.. -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums 
.. -fno-common -Wmissing-prototypes -Wnested-externs -Dinline= 
.. -g -O4 tmp.c 
.. cc1: warnings being treated as errors
.. tmp.c:7: warning: function declaration isn't a prototype
.. tmp.c: In function `main':
.. tmp.c:9: warning: passing arg 1 of `f' as integer rather than floating 
.. due to prototype
.. tmp.c: In function `f':
.. tmp.c:14: warning: type of `x' defaults to `int'
.. tmp.c:15: warning: passing arg 1 of `sqrt' as floating rather than integer 
.. due to prototype
.. make: *** [tmp] Error 1
.. @end example

.. @noindent
.. The error in the prototype is flagged, plus the fact that we should have
.. defined main as :code:`int main (void)` in ANSI C. Clearly there is some
.. work to do before this program is ready to run.

References and Further Reading
==============================

The following books are essential reading for anyone writing and
debugging numerical programs with :code:`gcc` and :code:`gdb`.

* R.M. Stallman, *Using and Porting GNU CC*, Free Software
  Foundation, ISBN 1882114388

* R.M. Stallman, R.H. Pesch, *Debugging with GDB: The GNU
  Source-Level Debugger*, Free Software Foundation, ISBN 1882114779

For a tutorial introduction to the GNU C Compiler and related programs,
see 

* B.J. Gough, http://www.network-theory.co.uk/gcc/intro/,'
  *An Introduction to GCC*, Network Theory
  Ltd, ISBN 0954161793
