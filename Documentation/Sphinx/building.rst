========
Building
========

.. contents:: Table of Contents
   :depth: 2
   :local:

Dependent Packages
==================

Required Packages
-----------------

-  `CMake <http://www.cmake.org>`__ (Version 2.8 or newer)
-  `Insight Toolkit`_ (ITK, version 4 or newer)
-  `Visualization Toolkit`_ (VTK, version 6 or newer)
-  `SlicerExecutionModel`_ (GenerateCLP, for command line interface):
-  `MKL`_ or `OpenBLAS`_ (at least one is required)
-  `GNU Scientific Library`_ (GSL, for special functions)

Optional Packages
-----------------

-  `Qt`_: needed for building ``QTApplications``.
-  `Doxygen`_: needed when building doxygen documents.

Third Party Packages (included in dmritool)
-------------------------------------------

-  `GTest`_: GTest can be automatically downloaded and built when building tests.
-  `SPAMS`_: included in ``ThirdParty/Spams``

Build Souce Code
================

Linux and Mac
--------------

1. build/install GSL_, VTK_, ITK_, Qt_ (optional), Doxygen_ (optional).
2. build/install MKL_ or OpenBlas_.

   MKL_ was free for non-commercial usage in linux. However now the 
   `free download link <https://software.intel.com/en-us/non-commercial-software-development>`__
   seems dead forever. 
   OpenBlas_ is always free and easy to use. 

   If you want to use OpenBlas_, please make sure that you build OpenBlas_ correctly
   or set related environment variables. 
   See `FAQ of OpenBlas <https://github.com/xianyi/OpenBLAS/wiki/faq#multi-threaded>`__.
   Please make sure set ``USE_THREAD=0 USE_OPENMP=1`` when building OpenBlas_. 
   Suggested building command for OpenBlas_ used in dmritool_ is: ::

        git clone https://github.com/xianyi/OpenBLAS.git
        cd OpenBlas
        make BINARY=64 USE_THREAD=0 USE_OPENMP=1 cc=gcc  FC=gfortran  
        make install

3. build SlicerExecutionModel_ ::

       git clone https://github.com/Slicer/SlicerExecutionModel.git
       mkdir SlicerExecutionModel-build
       cd SlicerExecutionModel-build
       ccmake ../SlicerExecutionModel
       make

4. build DMRITOOL_ source ::

       mkdir dmritool-build
       cd dmritool-build
       ccmake DIRECTORY_OF_DMRITOOL_SOURCE
       make

  -  When building DMRITOOL_ with ``ccmake``, ``GenerateCLP_DIR`` needs to be set as ``SlicerExecutionModel-build/GenerateCLP``.
  -  The code uses OpenBlas_ if ``DMRITOOL_USE_MKL`` is set off, it uses MKL_ if ``DMRITOOL_USE_MKL`` is ``ON``.
  -  If you want to build matlab_ mex files, please set ``DMRITOOL_WRAP_MATLAB`` ``ON``.
  -  GTest_ is needed when ``BUILD_TESTING`` is set ``ON`` for tests. We suggest you set it ``ON``
  -  Qt_ is needed if ``BUILD_QT_APPLICATIONS`` is set ``ON``.
  -  ``VERBOSITY_LEVEL`` is used for debug. The default value is 0. 
     If you set it as 1, the routines can provide more logging information and perform more condition checking, 
     but it may make the built executables a little bit slower.


5. build tests (optional, but suggested) ::

      make test

  - Set ``BUILD_TESTING`` ``ON`` to build tests. 
  - Make sure all tests are successfully passed. Otherwise the routines may not run correctly.  

6. build doxygen_ document (optional) ::

       make doxygen

7. Set environments. 

  - You can add ``dmritool-build/bin`` into your ``PATH`` to use the binary executables easily.
  - When set ``DMRITOOL_WRAP_MATLAB``, the built mex executables are in ``dmritool-build/Wrapping/Matlab/bin``. 
    Then you need to put the folder into matlab path. 
  - For matlab mex files, you also need to set ``BLAS_VERSION`` as the openblas lib or mkl lib. 
    If you use mkl_, then set ``blas_version`` as ::

       export BLAS_VERSION="MKL_DIR/lib/intel64/libmkl_rt.so"

    If you use openblas_, then set ``blas_version`` as ::

       export BLAS_VERSION="/usr/local/lib/libopenblas.so"

Windows
--------

We did not test the building of DMRITOOL_ in windows. 
You can try it by yourself. 



.. include:: links.inc
