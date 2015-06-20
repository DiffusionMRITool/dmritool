========================
DMRITOOL's documentation
========================

.. contents:: Table of Contents
   :depth: 2
   :local:

.. image:: https://travis-ci.org/DiffusionMRITool/dmritool.svg?branch=master
    :target: https://travis-ci.org/DiffusionMRITool/dmritool
    :alt: travis-ci Status             


Introduction
============

DMRITOOL_ is a **free** and **open source** toolbox for `diffusion MRI <http://en.wikipedia.org/wiki/Diffusion_MRI>`__ data processing. 
It is written in C++ with Matlab_ interface. 
Currently DMRITOOL_ has no GUI. It only provides executables for command line usage. 
You can also use the released mex executables in matlab_. 


With DMRITOOL_, you can: 

* perform reconstruction/estimation of diffusion data, including diffusion weighted signal, 
  ensemble average propagator (EAP), diffusion orientation distribution function (dODF), 
  and some meaningful scalar maps, etc.
* generate spherically uniform sampling schemes for single or multiple shells.
* perform diffusion MRI data simulation.  
* visualize spherical function field (e.g. dODF field, EAP profile field)


Citations will help us support the continued development of DMRITOOL_.
If you use the methods and codes released in DMRITOOL_, please cite the
related references. See :doc:`citation page <citation>`.

Download
========

You can download the latest source codes from `github <https://github.com/DiffusionMRITool/dmritool>`__::

  git clone https://github.com/DiffusionMRITool/dmritool

Acknowledgements
================

DMRITOOL_ is/was supported by the following research groups:

*  `Section on Tissue Biophysics and Biomimetics (STBB)`_, `National Institutes of Health`_
*  `MIND Lab`_, `IDEA Group`_, `University of North Carolina at Chapel Hill`_
*  `Athena project team`_, INRIA_
*  `Brainnetome center`_, `Institute of Automation, Chinese Academy of Sciences`_


Users Documentation
===================

.. toctree::
   :maxdepth: 2

   news
   building
   tutorials
   userguide
   commands/commandlist
   matlabfiles/matlabfunctions
   support
   citation
   license

Developers Documentation
========================

.. toctree::
   :maxdepth: 2

   developers
   doxygen
   codemodules


.. include:: links.inc
