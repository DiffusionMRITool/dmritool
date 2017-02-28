
:tocdepth: 1

====
News
====

.. meta::
   :description: news of dmritool
.. include:: meta_keywords.txt
.. include:: links.inc

.. contents:: Table of Contents
   :depth: 2
   :local:

2017-02-28: DMRITool now uses C++11
===================================

* We use C++11 in DMRITool now. 

  + NDArray now is more efficient with C++11. 
  + The codes have to be built with gcc (>= 4.8), ITK (>= 4.9). Please check :doc:`the building instruction <building>`. 

* Some algorithms and demos on sampling scheme design are released.

  + Please check the :doc:`updated tutorial on sampling scheme design<tutorial_qspacesampling>`.
  + 1-Opt: :doc:`SamplingSchemeQSpace1OptEstimation <commands/SamplingSchemeQSpace1OptEstimation>`, a new greedy algorithm which refines a given initialization. 
    It can fix a known issue of :doc:`IMOC <commands/SamplingSchemeQSpaceIMOCEstimation>` when the number of samples is large. 
    See `the tutorial <tutorial_qspacesampling.html#sampling-schemes-by-imoc-1-opt-cnlo>`__.
  + CNLO: :doc:`OptimalSamplingSingleShellCNLO <matlabfiles/OptimalSamplingSingleSubset>` for single-shell scheme design. 
  + CNLO: :doc:`OptimalSamplingMultiShellCNLO <matlabfiles/OptimalSamplingMultiShellCNLO>` for multiple-shell scheme design.
  + jointly use IMOC + 1-Opt + CNLO for sampling scheme design:
    
    + `demo for single shell scheme <demos/demo_singleshell_IMOC_1Opt_CNLO.html>`__
    + `demo for multiple shell scheme <demos/demo_mutishell_IMOC_1Opt_CNLO.html>`__



2015-06-18: Some codes on sampling scheme design are released
=============================================================

* Mixed Integer Linear Programming (MILP) for sampling scheme design is released in matlab. See folder ``${DMRITOOL_SOURCE_DIR}/Matlab``. 
* Iterative Maximum Overlap Construction (IMOC) for sampling scheme design is released in C++. 
  See command :doc:`SamplingSchemeQSpaceIMOCEstimation <commands/SamplingSchemeQSpaceIMOCEstimation>`.
* Add :doc:`a tutorial on sampling scheme design <tutorial_qspacesampling>` with some demos.
* Use some scripts to generate the c++ executable list, matlab function list, run demos in tutorials, update dmritool website and doxygen. 
* Add :doc:`VTKPolyData.py <commands/VTKPolyData.py>` to visualize vtk files and save the visualization into png files. It can replace :doc:`vtkviewer <commands/vtkviewer>`.


2015-06-09: DMRITool is using github pages
==========================================

* The homepage of DMRITool_ is moving from `readthedocs <http://dmritool.readthedocs.org/>`__ to github pages, 
* Using gh-pages, it is possible to automatically generate some documents, for example, :doc:`the list of commands <commands/commandlist>`. 
* As required by a user, :doc:`ComputeSHCoefficientsOfDWIFromSymmetricTensor <commands/ComputeSHCoefficientsOfDWIFromSymmetricTensor>` is added to compute a fiber response function from a given symmetric tensor. 
* Some tips using paraview_ is added based on the `discussion <http://www.nitrc.org/pipermail/dmritool-discussion/2015-March/000001.html>`__

2015-03-09: Version 0.1.1 is released
=====================================

The version v0.1.1 is released. 
It is mainly a bug-fix release to make the building process better. 

* We now use the dmritool-discussion_ mailing list by `nitrc <http://www.nitrc.org/projects/dmritool>`__. 
* The doxygen documentation can be found in the `github page <http://diffusionmritool.github.io/dmritool-doxygen>`__. 

2015-02-18: DMRITool is now open source
=======================================

DMRITool_ is released as an open source toolbox. 
The source code is in `github <https://github.com/DiffusionMRITool/dmritool>`__ and the documentation is in `readthedocs <http://dmritool.readthedocs.org/>`__.

The fist release (v0.1) includes:

* Spherical Polar Fourier Imaging. See :doc:`the SPFI tutorial <tutorial_spfi>`.
* DWI simulation. See :doc:`the DWI data simulation tutotial <tutorial_dwisimulation>`.
* Visualization of spherical function fields.



