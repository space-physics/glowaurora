================
Quickstart Guide
================

1. Print out files QUICKSTART, README and glow.txt
2. Copy all the files to a directory
3. Copy nrlmsise00.f and iri90.f to the directory or tell the compiler where they are
4. Copy the IRI data files to another directory
5. Rename the example program dayexample.driver to dayexample.f
6. Edit the IRI call in dayexample.f to specify the IRI data file directory
7. Compile and link *.f in Fortran 77 with -O3 optimization

Daytime Example
----------------
1. Run the executable (a.out) with input from dayexample.in, i.e.::

    a.out < dayexample.in > test1.out

2. Compare the output you get to dayexample.out.  If it looks reasonable, try making some plots

Aurora Example
--------------
1. For auroral runs, change aurexample.driver to aurexample.f, compile, link::

    a.out < aurexample.in > test2.out

2. Compare the output you get to aurexample.out

High-energy aurora example
--------------------------
1. change the name of hexexample.driver to hexexample.f
2. edit the glow.h file:
    a)  comment out the standard parameters
    b)  comment in the high-energy parameters

3. re-compile all subroutines, link, and::

    a.out < hexexample.in > test3.out

4. compare to hexexample.out

After that
----------
Try some runs with your own input parameters.

Then try making some modifications to the example program to suit your purposes

You can also run it all out of IDL using:

.. code-block:: idl

  spawn, './a.out < input.file > output.file'
