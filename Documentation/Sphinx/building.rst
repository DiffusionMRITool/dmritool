========
Building
========

.. meta::
   :description: building of dmritool
.. include:: meta_keywords.txt
.. include:: links.inc

.. contents:: Table of Contents
   :depth: 2
   :local:

.. image:: https://travis-ci.org/DiffusionMRITool/dmritool.svg?branch=master
    :target: https://travis-ci.org/DiffusionMRITool/dmritool

Matlab codes
============

The codes in folder ``${DMRITOOL_SOURCE_DIR}/Matlab`` are purely in matlab. 
Thus you do not need to build c++ codes to run codes in that folder. 
The codes in folder ``${DMRITOOL_SOURCE_DIR}/Wrapping/Matlab`` are matlab mex files, which depend on the c++ codes. 

See :doc:`the matlab function list <matlabfiles/matlabfunctions>`, 
where these functions starting with mex prefix need to be built using cmake. 

Dependent Packages
==================

Required Packages
-----------------

-  CMake_ (Version 2.8 or newer)
-  GCC_ (4.8 or newer. DMRITool uses C++11, thus it cannot work for old gcc)
-  `Insight Toolkit`_ (ITK, version 4.9 or newer)
-  `Visualization Toolkit`_ (VTK, version 6 or newer)
-  SlicerExecutionModel_ (GenerateCLP, for command line interface):
-  `OpenBLAS`_ + Lapack_ or MKL_ (at least one is required)
-  `GNU Scientific Library`_ (GSL, for special functions)

Optional Packages
-----------------

-  Qt_: needed for building ``QTApplications``.
-  Doxygen_: needed when building doxygen documents.

Third Party Packages (included in dmritool)
-------------------------------------------

-  GTest_: GTest can be automatically downloaded and built when building tests.
-  SPAMS_: included in ``ThirdParty/Spams``

Build Souce Code
================

Linux and Mac
--------------

1. build/install CMake_, GCC_, Git_, GSL_, VTK_, ITK_, Qt_ (optional), Doxygen_ (optional).
2. build/install OpenBlas_ + Lapack_ or MKL_ .

   MKL_ is free for students, educators, academic researchers and open source contributors. 
   See `the free download link <https://software.intel.com/en-us/qualify-for-free-software>`__. 
   OpenBlas_ and Lapack_ are always free and easy to use. 

   If you choose OpenBlas_ + Lapack_, you can use them from system repositories. 
   In this case, you need to set environment variables correctly. 
   See `FAQ of OpenBlas <https://github.com/xianyi/OpenBLAS/wiki/faq#multi-threaded>`__.
   ::

        export OPENBLAS_NUM_THREADS=1
    
   We suggest that you build OpenBlas_ from its source codes. 
   Suggested building command for OpenBlas_ used in dmritool_ is: ::

        git clone https://github.com/xianyi/OpenBLAS.git
        cd OpenBlas
        git checkout v0.2.13
        make BINARY=64 USE_THREAD=0 USE_OPENMP=1 cc=gcc  FC=gfortran  
        make install

3. build SlicerExecutionModel_ ::

       git clone https://github.com/Slicer/SlicerExecutionModel.git
       mkdir SlicerExecutionModel-build
       cd SlicerExecutionModel-build
       ccmake ../SlicerExecutionModel
       make

4. build DMRITool_ source ::

       git clone https://github.com/DiffusionMRITool/dmritool
       mkdir dmritool-build
       cd dmritool-build
       ccmake ../dmritool
       make

  -  When building DMRITool_ with ``ccmake``, ``GenerateCLP_DIR`` needs to be set as ``SlicerExecutionModel-build/GenerateCLP``.
  -  The code uses OpenBlas_ + Lapack_ if ``DMRITOOL_USE_MKL`` is set ``OFF``, it uses MKL_ if ``DMRITOOL_USE_MKL`` is ``ON``.
  -  If ``DMRITOOL_USE_FASTLAPACK=ON``, then the codes use fast versions for SVD and eigen-decomposition, 
     but it may have errors for openblas_. 
     You can check whether ``utlVNLLapackGTest`` and ``utlVNLBlasGTest`` can successfully pass. 
     If these two tests do not passed, you may need to build OpenBlas_ manually or build DMRITool_ with ``DMRITOOL_USE_FASTLAPACK=OFF``
  -  If you want to build matlab_ mex files, please set ``DMRITOOL_WRAP_MATLAB=ON``.
  -  If ``BUILD_TESTING`` is set ``ON``, some tests based on GTest_ will be built. 
     We suggest you set it ``ON``. GTest_ will be automatically downloaded and built if you do not have it in system. 
  -  Qt_ is needed if ``BUILD_QT_APPLICATIONS`` is set ``ON``.
  -  ``VERBOSITY_LEVEL`` is used for debug. The default value is 0. 
     If you set it as 1, the routines can provide more logging information and perform more condition checking, 
     but it may make the built executables a little bit slower.
  -  If you use an old version of ITK (<4.9), you may have to set ``CMAKE_CXX_FLAGS=-fpermissive`` to build the codes. 
     That is `a known issue of ITK with c++11 support <https://itk.org/pipermail/insight-developers/2015-December/024705.html>`_, which is solved after ITK 4.9. 

5. build doxygen_ document (optional) ::

      make doxygen

6. Set environments. 

  - You can add ``dmritool-build/bin`` into your ``PATH`` to use the binary executables easily.
  - When ``DMRITOOL_WRAP_MATLAB`` is set ``ON``, the built mex executables are in ``dmritool-build/Wrapping/Matlab/bin``. 
    Then to use these mex executables in matlab, you need to put the folder into matlab path. 
  - For matlab mex files, you also need to set ``BLAS_VERSION`` as the openblas lib or mkl lib. 
    If you use mkl_, then set ``blas_version`` as ::

       export BLAS_VERSION="MKL_DIR/lib/intel64/libmkl_rt.so"

    If you use openblas_, then set ``blas_version`` as ::

       export BLAS_VERSION="/usr/local/lib/libopenblas.so"

7. build tests (optional, but suggested) ::

      make test

  - Set ``BUILD_TESTING`` ``ON`` to build tests. 
  - Make sure all tests are successfully passed. Otherwise the routines may not run correctly.  


Windows
--------

We did not test the building of DMRITool_ in windows. 
You can try it by yourself. 



