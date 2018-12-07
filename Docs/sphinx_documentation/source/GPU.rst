.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


.. _sec:gpu:build:

Building GPU Support
====================

Building with GNU Make
----------------------

To build AMReX with GPU support, you add ``USE_CUDA=TRUE`` to your
make file or as an command line argument.  AMReX itself does not
require OpenACC and CUDA Fortran, but application codes can use them
if they are supported by the compiler.  For OpenACC support, you add
``USE_ACC=TRUE``.  Only IBM and PGI support CUDA Fortran.  Thus for
CUDA Fortran, you must use ``COMP=pgi`` or ``COMP=ibm``.  The default
host compiler for NVCC is GCC even if ``COMP`` is set to a different
compiler.  One can change this by setting ``NVCC_HOST_COMP`` to say
``pgi``.  For example, ``COMP=pgi`` alone will compile C/C++ codes
with NVCC/GCC and Fortran codes with PGI, and link with PGI.  Using
``COMP=pgi`` and ``NVCC_HOST_COMP=pgi`` will compile C/C++ codes with
NVCC/PGI.

You can use ``Tutorials/Basic/HelloWorld_C`` to test your programming
environment.  Typing

.. highlight:: console

::

   make COMP=gnu USE_CUDA=TRUE

should produce an executable named ``main3d.gnu.DEBUG.CUDA.ex``.  You
can run it and that will generate results like below

.. highlight:: console

::

   $ ./main3d.gnu.DEBUG.CUDA.ex 
   CUDA initialized with 1 GPU
   AMReX (18.12-95-gf265b537f479-dirty) initialized
   Hello world from AMReX version 18.12-95-gf265b537f479-dirty
   [The         Arena] space (kilobyte): 8192
   [The  Device Arena] space (kilobyte): 8192
   [The Managed Arena] space (kilobyte): 8192
   [The  Pinned Arena] space (kilobyte): 8192
   AMReX (18.12-95-gf265b537f479-dirty) finalized

.. _sec:gpu:namespace:

Gpu Namespace and Macros
========================


.. _sec:gpu:memory:

Memory Allocation
=================

memory


.. _sec:gpu:launch:

Kernel Launch
=============

how to launch a kernel

streams and mfiter

inLaunchRegion

macro, cuda fortran and openacc

refer back to Basic for calling C++ functions on FArrayBox


.. _sec:gpu:reduction:

Reduction
=========


.. _sec:gpu:asyncfab:

AsyncFab and AsyncArray
=======================


.. _sec:gpu:particle:

Particle
========


.. _sec::gpu:mpi:

CUDA Aware MPI
==============


.. _sec:gpu:limits:

Limitations
===========

At most one gpu per mpi rank.  OpenMP, AMR development are underway
cmake and build as library
