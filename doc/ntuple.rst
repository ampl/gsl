.. index:: ntuples

********
N-tuples
********

This chapter describes functions for creating and manipulating
*ntuples*, sets of values associated with events.  The ntuples
are stored in files. Their values can be extracted in any combination
and *booked* in a histogram using a selection function.

The values to be stored are held in a user-defined data structure, and
an ntuple is created associating this data structure with a file.  The
values are then written to the file (normally inside a loop) using
the ntuple functions described below.

A histogram can be created from ntuple data by providing a selection
function and a value function.  The selection function specifies whether
an event should be included in the subset to be analyzed or not. The value
function computes the entry to be added to the histogram for each
event.

All the ntuple functions are defined in the header file
:file:`gsl_ntuple.h`.

The ntuple struct
=================

.. type:: gsl_ntuple

   Ntuples are manipulated using the :type:`gsl_ntuple` struct. This struct
   contains information on the file where the ntuple data is stored, a
   pointer to the current ntuple data row and the size of the user-defined
   ntuple data struct::

      typedef struct
        {
          FILE * file;
          void * ntuple_data;
          size_t size;
        } gsl_ntuple;

Creating ntuples
================

.. function:: gsl_ntuple * gsl_ntuple_create (char * filename, void * ntuple_data, size_t size)

   This function creates a new write-only ntuple file :data:`filename` for
   ntuples of size :data:`size` and returns a pointer to the newly created
   ntuple struct.  Any existing file with the same name is truncated to
   zero length and overwritten.  A pointer to memory for the current ntuple
   row :data:`ntuple_data` must be supplied---this is used to copy ntuples
   in and out of the file.

Opening an existing ntuple file
===============================

.. function:: gsl_ntuple * gsl_ntuple_open (char * filename, void * ntuple_data, size_t size)

   This function opens an existing ntuple file :data:`filename` for reading
   and returns a pointer to a corresponding ntuple struct. The ntuples in
   the file must have size :data:`size`.  A pointer to memory for the current
   ntuple row :data:`ntuple_data` must be supplied---this is used to copy
   ntuples in and out of the file.

Writing ntuples
===============

.. function:: int gsl_ntuple_write (gsl_ntuple * ntuple)

   This function writes the current ntuple :code:`ntuple->ntuple_data` of
   size :code:`ntuple->size` to the corresponding file.

.. function:: int gsl_ntuple_bookdata (gsl_ntuple * ntuple)

   This function is a synonym for :func:`gsl_ntuple_write`.

Reading ntuples
===============

.. function:: int gsl_ntuple_read (gsl_ntuple * ntuple)

   This function reads the current row of the ntuple file for :data:`ntuple`
   and stores the values in :code:`ntuple->data`.

Closing an ntuple file
======================

.. function:: int gsl_ntuple_close (gsl_ntuple * ntuple)

   This function closes the ntuple file :data:`ntuple` and frees its
   associated allocated memory.

Histogramming ntuple values
===========================

Once an ntuple has been created its contents can be histogrammed in
various ways using the function :func:`gsl_ntuple_project`.  Two
user-defined functions must be provided, a function to select events and
a function to compute scalar values. The selection function and the
value function both accept the ntuple row as a first argument and other
parameters as a second argument.

.. index::
   single: selection function, ntuples

.. type:: gsl_ntuple_select_fn

   The *selection function* determines which ntuple rows are selected
   for histogramming.  It is defined by the following struct::

      typedef struct
        {
          int (* function) (void * ntuple_data, void * params);
          void * params;
        } gsl_ntuple_select_fn;

   The struct component :data:`function` should return a non-zero value for
   each ntuple row that is to be included in the histogram.

.. index::
   single: value function, ntuples

.. type:: gsl_ntuple_value_fn

   The *value function* computes scalar values for those ntuple rows
   selected by the selection function::

      typedef struct
        {
          double (* function) (void * ntuple_data, void * params);
          void * params;
        } gsl_ntuple_value_fn;

   In this case the struct component :data:`function` should return the value
   to be added to the histogram for the ntuple row.  

.. index::
   single: histogram, from ntuple
   single: projection of ntuples

.. function:: int gsl_ntuple_project (gsl_histogram * h, gsl_ntuple * ntuple, gsl_ntuple_value_fn * value_func, gsl_ntuple_select_fn * select_func)

   This function updates the histogram :data:`h` from the ntuple :data:`ntuple`
   using the functions :data:`value_func` and :data:`select_func`. For each
   ntuple row where the selection function :data:`select_func` is non-zero the
   corresponding value of that row is computed using the function
   :data:`value_func` and added to the histogram.  Those ntuple rows where
   :data:`select_func` returns zero are ignored.  New entries are added to
   the histogram, so subsequent calls can be used to accumulate further
   data in the same histogram.

Examples
========

The following example programs demonstrate the use of ntuples in
managing a large dataset.  The first program creates a set of 10,000
simulated "events", each with 3 associated values :math:`(x,y,z)`.  These
are generated from a Gaussian distribution with unit variance, for
demonstration purposes, and written to the ntuple file :file:`test.dat`.

.. include:: examples/ntuplew.c
   :code:

The next program analyses the ntuple data in the file :file:`test.dat`.
The analysis procedure is to compute the squared-magnitude of each
event, :math:`E^2=x^2+y^2+z^2`, and select only those which exceed a
lower limit of 1.5.  The selected events are then histogrammed using
their :math:`E^2` values.

.. include:: examples/ntupler.c
   :code:

:numref:`fig_ntuples` shows the distribution of the selected events.
Note the cut-off at the lower bound.

.. _fig_ntuples:

.. figure:: /images/ntuple.png
   :scale: 60%

   Distribution of selected events

References and Further Reading
==============================

.. index:: PAW, HBOOK

Further information on the use of ntuples can be found in the
documentation for the CERN packages PAW and HBOOK
(available online).
