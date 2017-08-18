.. index::
   single: physical constants
   single: constants, physical
   single: conversion of units
   single: units, conversion of

******************
Physical Constants
******************

This chapter describes macros for the values of physical constants, such
as the speed of light, :math:`c`, and gravitational constant, :math:`G`.
The values are available in different unit systems, including the
standard MKSA system (meters, kilograms, seconds, amperes) and the CGSM
system (centimeters, grams, seconds, gauss), which is commonly used in
Astronomy.

The definitions of constants in the MKSA system are available in the file
:file:`gsl_const_mksa.h`.  The constants in the CGSM system are defined in
:file:`gsl_const_cgsm.h`.  Dimensionless constants, such as the fine
structure constant, which are pure numbers are defined in
:file:`gsl_const_num.h`.

The full list of constants is described briefly below.  Consult the
header files themselves for the values of the constants used in the
library.

.. index::
   single: fundamental constants
   single: constants, fundamental

Fundamental Constants
=====================

.. macro:: GSL_CONST_MKSA_SPEED_OF_LIGHT

   The speed of light in vacuum, :math:`c`.

.. macro:: GSL_CONST_MKSA_VACUUM_PERMEABILITY

   The permeability of free space, :math:`\mu_0`. This constant is defined
   in the MKSA system only.

.. macro:: GSL_CONST_MKSA_VACUUM_PERMITTIVITY

   The permittivity of free space, :math:`\epsilon_0`.  This constant is
   defined in the MKSA system only.

.. macro:: GSL_CONST_MKSA_PLANCKS_CONSTANT_H

   Planck's constant, :math:`h`.

.. macro:: GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR

   Planck's constant divided by :math:`2\pi`, :math:`\hbar`.

.. macro:: GSL_CONST_NUM_AVOGADRO

   Avogadro's number, :math:`N_a`.

.. macro:: GSL_CONST_MKSA_FARADAY

   The molar charge of 1 Faraday.

.. macro:: GSL_CONST_MKSA_BOLTZMANN

   The Boltzmann constant, :math:`k`.

.. macro:: GSL_CONST_MKSA_MOLAR_GAS

   The molar gas constant, :math:`R_0`.

.. macro:: GSL_CONST_MKSA_STANDARD_GAS_VOLUME

   The standard gas volume, :math:`V_0`.

.. macro:: GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT

   The Stefan-Boltzmann radiation constant, :math:`\sigma`.

.. macro:: GSL_CONST_MKSA_GAUSS

   The magnetic field of 1 Gauss.

.. index:: astronomical constants

Astronomy and Astrophysics
==========================

.. macro:: GSL_CONST_MKSA_ASTRONOMICAL_UNIT

   The length of 1 astronomical unit (mean earth-sun distance), :math:`au`.

.. macro:: GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT

   The gravitational constant, :math:`G`.

.. macro:: GSL_CONST_MKSA_LIGHT_YEAR

   The distance of 1 light-year, :math:`ly`.

.. macro:: GSL_CONST_MKSA_PARSEC

   The distance of 1 parsec, :math:`pc`.

.. macro:: GSL_CONST_MKSA_GRAV_ACCEL

   The standard gravitational acceleration on Earth, :math:`g`.

.. macro:: GSL_CONST_MKSA_SOLAR_MASS

   The mass of the Sun.

.. index::
   single: atomic physics, constants
   single: nuclear physics, constants

Atomic and Nuclear Physics
==========================

.. macro:: GSL_CONST_MKSA_ELECTRON_CHARGE

   The charge of the electron, :math:`e`.

.. macro:: GSL_CONST_MKSA_ELECTRON_VOLT

   The energy of 1 electron volt, :math:`eV`.

.. macro:: GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS

   The unified atomic mass, :math:`amu`.

.. macro:: GSL_CONST_MKSA_MASS_ELECTRON

   The mass of the electron, :math:`m_e`.

.. macro:: GSL_CONST_MKSA_MASS_MUON

   The mass of the muon, :math:`m_\mu`.

.. macro:: GSL_CONST_MKSA_MASS_PROTON

   The mass of the proton, :math:`m_p`.

.. macro:: GSL_CONST_MKSA_MASS_NEUTRON

   The mass of the neutron, :math:`m_n`.

.. macro:: GSL_CONST_NUM_FINE_STRUCTURE

   The electromagnetic fine structure constant :math:`\alpha`.

.. macro:: GSL_CONST_MKSA_RYDBERG

   The Rydberg constant, :math:`Ry`, in units of energy.  This is related to
   the Rydberg inverse wavelength :math:`R_\infty` by :math:`Ry = h c R_\infty`.

.. macro:: GSL_CONST_MKSA_BOHR_RADIUS

   The Bohr radius, :math:`a_0`.

.. macro:: GSL_CONST_MKSA_ANGSTROM

   The length of 1 angstrom.

.. macro:: GSL_CONST_MKSA_BARN

   The area of 1 barn.

.. macro:: GSL_CONST_MKSA_BOHR_MAGNETON

   The Bohr Magneton, :math:`\mu_B`.

.. macro:: GSL_CONST_MKSA_NUCLEAR_MAGNETON

   The Nuclear Magneton, :math:`\mu_N`.

.. macro:: GSL_CONST_MKSA_ELECTRON_MAGNETIC_MOMENT

   The absolute value of the magnetic moment of the electron, :math:`\mu_e`.
   The physical magnetic moment of the electron is negative.

.. macro:: GSL_CONST_MKSA_PROTON_MAGNETIC_MOMENT

   The magnetic moment of the proton, :math:`\mu_p`.

.. macro:: GSL_CONST_MKSA_THOMSON_CROSS_SECTION

   The Thomson cross section, :math:`\sigma_T`.

.. macro:: GSL_CONST_MKSA_DEBYE

   The electric dipole moment of 1 Debye, :math:`D`.

.. index:: time units

Measurement of Time
===================

.. macro:: GSL_CONST_MKSA_MINUTE

   The number of seconds in 1 minute.

.. macro:: GSL_CONST_MKSA_HOUR

   The number of seconds in 1 hour.

.. macro:: GSL_CONST_MKSA_DAY

   The number of seconds in 1 day.

.. macro:: GSL_CONST_MKSA_WEEK

   The number of seconds in 1 week.

.. index::
   single: imperial units
   single: units, imperial

Imperial Units
==============

.. macro:: GSL_CONST_MKSA_INCH

   The length of 1 inch.

.. macro:: GSL_CONST_MKSA_FOOT

   The length of 1 foot.

.. macro:: GSL_CONST_MKSA_YARD

   The length of 1 yard.

.. macro:: GSL_CONST_MKSA_MILE

   The length of 1 mile.

.. macro:: GSL_CONST_MKSA_MIL

   The length of 1 mil (1/1000th of an inch).

.. index:: nautical units

Speed and Nautical Units
========================

.. macro:: GSL_CONST_MKSA_KILOMETERS_PER_HOUR

   The speed of 1 kilometer per hour.

.. macro:: GSL_CONST_MKSA_MILES_PER_HOUR

   The speed of 1 mile per hour.

.. macro:: GSL_CONST_MKSA_NAUTICAL_MILE

   The length of 1 nautical mile.

.. macro:: GSL_CONST_MKSA_FATHOM

   The length of 1 fathom.

.. macro:: GSL_CONST_MKSA_KNOT

   The speed of 1 knot.

.. index:: printers units

Printers Units
==============

.. macro:: GSL_CONST_MKSA_POINT

   The length of 1 printer's point (1/72 inch).

.. macro:: GSL_CONST_MKSA_TEXPOINT

   The length of 1 TeX point (1/72.27 inch).

.. index:: volume units

Volume, Area and Length
=======================

.. macro:: GSL_CONST_MKSA_MICRON

   The length of 1 micron.

.. macro:: GSL_CONST_MKSA_HECTARE

   The area of 1 hectare.

.. macro:: GSL_CONST_MKSA_ACRE

   The area of 1 acre.

.. macro:: GSL_CONST_MKSA_LITER

   The volume of 1 liter.

.. macro:: GSL_CONST_MKSA_US_GALLON

   The volume of 1 US gallon.

.. macro:: GSL_CONST_MKSA_CANADIAN_GALLON

   The volume of 1 Canadian gallon.

.. macro:: GSL_CONST_MKSA_UK_GALLON

   The volume of 1 UK gallon.

.. macro:: GSL_CONST_MKSA_QUART

   The volume of 1 quart.

.. macro:: GSL_CONST_MKSA_PINT

   The volume of 1 pint.

.. @node Cookery
.. @section Cookery
.. @commentindex cookery units

.. @table @commentode
.. @item GSL_CONST_MKSA_CUP
.. The volume of 1 cup.

.. @item GSL_CONST_MKSA_FLUID_OUNCE
.. The volume of 1 fluid ounce.

.. @item GSL_CONST_MKSA_TABLESPOON
.. The volume of 1 tablespoon.

.. @item GSL_CONST_MKSA_TEASPOON
.. The volume of 1 teaspoon.
.. @end table

.. index::
   single: mass, units of
   single: weight, units of

Mass and Weight
===============

.. macro:: GSL_CONST_MKSA_POUND_MASS

   The mass of 1 pound.

.. macro:: GSL_CONST_MKSA_OUNCE_MASS

   The mass of 1 ounce.

.. macro:: GSL_CONST_MKSA_TON

   The mass of 1 ton.

.. macro:: GSL_CONST_MKSA_METRIC_TON

   The mass of 1 metric ton (1000 kg).

.. macro:: GSL_CONST_MKSA_UK_TON

   The mass of 1 UK ton.

.. macro:: GSL_CONST_MKSA_TROY_OUNCE

   The mass of 1 troy ounce.

.. macro:: GSL_CONST_MKSA_CARAT

   The mass of 1 carat.

.. macro:: GSL_CONST_MKSA_GRAM_FORCE

   The force of 1 gram weight.

.. macro:: GSL_CONST_MKSA_POUND_FORCE

   The force of 1 pound weight.

.. macro:: GSL_CONST_MKSA_KILOPOUND_FORCE

   The force of 1 kilopound weight.

.. macro:: GSL_CONST_MKSA_POUNDAL

   The force of 1 poundal.

.. index::
   single: energy, units of
   single: power, units of
   single: thermal energy, units of

Thermal Energy and Power
========================

.. macro:: GSL_CONST_MKSA_CALORIE

   The energy of 1 calorie.

.. macro:: GSL_CONST_MKSA_BTU

   The energy of 1 British Thermal Unit, :math:`btu`.

.. macro:: GSL_CONST_MKSA_THERM

   The energy of 1 Therm.

.. macro:: GSL_CONST_MKSA_HORSEPOWER

   The power of 1 horsepower.

.. index:: pressure, units of

Pressure
========

.. macro:: GSL_CONST_MKSA_BAR

   The pressure of 1 bar.

.. macro:: GSL_CONST_MKSA_STD_ATMOSPHERE

   The pressure of 1 standard atmosphere.

.. macro:: GSL_CONST_MKSA_TORR

   The pressure of 1 torr.

.. macro:: GSL_CONST_MKSA_METER_OF_MERCURY

   The pressure of 1 meter of mercury.

.. macro:: GSL_CONST_MKSA_INCH_OF_MERCURY

   The pressure of 1 inch of mercury.

.. macro:: GSL_CONST_MKSA_INCH_OF_WATER

   The pressure of 1 inch of water.

.. macro:: GSL_CONST_MKSA_PSI

   The pressure of 1 pound per square inch.

.. index:: viscosity, units of

Viscosity
=========

.. macro:: GSL_CONST_MKSA_POISE

   The dynamic viscosity of 1 poise.

.. macro:: GSL_CONST_MKSA_STOKES

   The kinematic viscosity of 1 stokes.

.. index::
   single: light, units of
   single: illumination, units of

Light and Illumination
======================

.. macro:: GSL_CONST_MKSA_STILB

   The luminance of 1 stilb.

.. macro:: GSL_CONST_MKSA_LUMEN

   The luminous flux of 1 lumen.

.. macro:: GSL_CONST_MKSA_LUX

   The illuminance of 1 lux.

.. macro:: GSL_CONST_MKSA_PHOT

   The illuminance of 1 phot.

.. macro:: GSL_CONST_MKSA_FOOTCANDLE

   The illuminance of 1 footcandle.

.. macro:: GSL_CONST_MKSA_LAMBERT

   The luminance of 1 lambert.

.. macro:: GSL_CONST_MKSA_FOOTLAMBERT

   The luminance of 1 footlambert.

.. index:: radioactivity, units of

Radioactivity
=============

.. macro:: GSL_CONST_MKSA_CURIE

   The activity of 1 curie.

.. macro:: GSL_CONST_MKSA_ROENTGEN

   The exposure of 1 roentgen.

.. macro:: GSL_CONST_MKSA_RAD

   The absorbed dose of 1 rad.

.. index:: force and energy, units of

Force and Energy
================

.. macro:: GSL_CONST_MKSA_NEWTON

   The SI unit of force, 1 Newton.

.. macro:: GSL_CONST_MKSA_DYNE

   The force of 1 Dyne = :math:`10^{-5}` Newton.

.. macro:: GSL_CONST_MKSA_JOULE

   The SI unit of energy, 1 Joule.

.. macro:: GSL_CONST_MKSA_ERG 

   The energy 1 erg = :math:`10^{-7}` Joule.

.. index::
   single: prefixes
   single: constants, prefixes

Prefixes
========

These constants are dimensionless scaling factors.

.. macro:: GSL_CONST_NUM_YOTTA

   :math:`10^{24}`

.. macro:: GSL_CONST_NUM_ZETTA

   :math:`10^{21}`

.. macro:: GSL_CONST_NUM_EXA

   :math:`10^{18}`

.. macro:: GSL_CONST_NUM_PETA

   :math:`10^{15}`

.. macro:: GSL_CONST_NUM_TERA

   :math:`10^{12}`

.. macro:: GSL_CONST_NUM_GIGA

   :math:`10^9`

.. macro:: GSL_CONST_NUM_MEGA

   :math:`10^6`

.. macro:: GSL_CONST_NUM_KILO

   :math:`10^3`

.. macro:: GSL_CONST_NUM_MILLI

   :math:`10^{-3}`

.. macro:: GSL_CONST_NUM_MICRO

   :math:`10^{-6}`

.. macro:: GSL_CONST_NUM_NANO

   :math:`10^{-9}`

.. macro:: GSL_CONST_NUM_PICO

   :math:`10^{-12}`

.. macro:: GSL_CONST_NUM_FEMTO

   :math:`10^{-15}`

.. macro:: GSL_CONST_NUM_ATTO
 
   :math:`10^{-18}`

.. macro:: GSL_CONST_NUM_ZEPTO

   :math:`10^{-21}`

.. macro:: GSL_CONST_NUM_YOCTO

   :math:`10^{-24}`

Examples
========

The following program demonstrates the use of the physical constants in
a calculation.  In this case, the goal is to calculate the range of
light-travel times from Earth to Mars.

The required data is the average distance of each planet from the Sun in
astronomical units (the eccentricities and inclinations of the orbits
will be neglected for the purposes of this calculation).  The average
radius of the orbit of Mars is 1.52 astronomical units, and for the
orbit of Earth it is 1 astronomical unit (by definition).  These values
are combined with the MKSA values of the constants for the speed of
light and the length of an astronomical unit to produce a result for the
shortest and longest light-travel times in seconds.  The figures are
converted into minutes before being displayed.

.. include:: examples/const.c
   :code:

Here is the output from the program,

.. include:: examples/const.txt
   :code:

References and Further Reading
==============================

The authoritative sources for physical constants are the 2006 CODATA
recommended values, published in the article below. Further
information on the values of physical constants is also available from
the NIST website.

* P.J. Mohr, B.N. Taylor, D.B. Newell, "CODATA Recommended
  Values of the Fundamental Physical Constants: 2006", Reviews of
  Modern Physics, 80(2), pp. 633--730 (2008).

* http://www.physics.nist.gov/cuu/Constants/index.html

* http://physics.nist.gov/Pubs/SP811/appenB9.html
