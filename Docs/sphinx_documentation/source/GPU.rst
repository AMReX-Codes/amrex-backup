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
CUDA Fortran, you must use ``COMP=pgi`` or ``COMP=ibm``.  OpenMP is
currently not supported when ``USE_CUDA=TRUE``.  The default host
compiler for NVCC is GCC even if ``COMP`` is set to a different
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

.. ===================================================================

.. _sec:gpu:namespace:

Gpu Namespace and Macros
========================

GPU related classes and functions are usually in ``namespace Gpu``,
which is inside ``namespace amrex``.  For portability, AMReX defines
some macros for CUDA function qualifiers and they should be preferred
to hardwired ``__*__``.  These include

.. highlight:: c++

::

   #define AMREX_GPU_HOST        __host__
   #define AMREX_GPU_DEVICE      __device__
   #define AMREX_GPU_GLOBAL      __global__
   #define AMREX_GPU_HOST_DEVICE __host__ __device__

Note that when AMReX is not built with CUDA, these macros expand to
empty space.

When AMReX is compiled with ``USE_CUDA=TRUE``, we pass
``-DAMREX_USE_CUDA`` and ``-DAMREX_USE_GPU`` to the compiler so that
these macros can be used for preprocessing.  For PGI and IBM
compilers, we also pass ``-DAMREX_USE_CUDA_FORTRAN``,
``-DAMREX_CUDA_FORT_GLOBAL='attributes(global)'``,
``-DAMREX_CUDA_FORT_DEVICE='attributes(device)'``, and
``-DAMREX_CUDA_FORT_HOST='attributes(host)'`` so that CUDA Fortran
functions can be properly declared.  When AMReX is compiled with
``USE_ACC=TRUE``, we pass ``-DAMREX_USE_ACC`` to the compiler.

.. ===================================================================

.. _sec:gpu:memory:

Memory Allocation
=================

To provide portability and improve memory allocation performance,
AMReX provides a number of memory pools.

.. raw:: latex

    \begin{center}

.. _tab:gpu:arena:

.. table:: Memory Arenas

    +---------------------+------------------+
    | Arena               |    Memory Type   |
    +=====================+==================+
    | The_Arena()         |  unified memory  | 
    +---------------------+------------------+
    | The_Device_Arena()  |  device memory   | 
    +---------------------+------------------+
    | The_Managed_Arena() |  unified memory  | 
    +---------------------+------------------+
    | The_Pinned_Arena()  |  pinned memory   | 
    +---------------------+------------------+

.. raw:: latex

    \end{center}

``Arena`` object returned by these arena functions provides

.. highlight:: c++

::

   void* alloc (std::size_t sz);
   void free (void* p);

``The_Arena()`` is used for memory allocation of data in ``BaseFab``.
Therefore the data in a ``MultiFab`` are in unified memory.
``The_Managed_Arena()`` is also a memory pool of unified memory, but
it is separated from ``The_Arena()`` for performance reason.  If you
want to print out the current memory usage by the arenas, you can call
``amrex::Arena::PrintUsage()``.

.. ===================================================================

.. _sec:gpu:launch:

Kernel Launch
=============

AMReX uses ``MFIter`` to iterate over a ``MultiFab``.  Inside the
loop, we call functions to work on ``FArrayBox`` objects
(see :ref:`sec:basics:mfiter`).  With GPU, we launch kernels inside
``MFIter`` loop.  A tutorial example can be found in
``Tutorials/GPU/Launch``.  The kernel launch part is shown below.



streams and mfiter

inLaunchRegion

macro, cuda fortran and openacc

refer back to Basic for calling C++ functions on FArrayBox

.. ===================================================================

.. _sec:gpu:safeclasses:

GPU Safe Classes
================

.. ===================================================================

.. _sec:gpu:assertion:

Assertion and Error Check
=========================

AMREX_GPU_SAFE_CALL(), AMREX_GPU_ERROR_CHECK();

.. ===================================================================

.. _sec:gpu:reduction:

Reduction
=========

.. ===================================================================

.. _sec:gpu:particle:

Particle
========

.. ===================================================================

.. _sec::gpu:mpi:

CUDA Aware MPI
==============

.. ===================================================================

.. _sec:gpu:limits:

Limitations
===========

At most one gpu per mpi rank.  OpenMP, AMR development are underway
cmake and build as library
