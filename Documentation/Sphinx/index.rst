====================================
DMRITOOL's documentation
====================================


Introduction
============

DMRITOOL_ is a **free** and **open source** toolbox for `diffusion
MRI <http://en.wikipedia.org/wiki/Diffusion_MRI>`__ data processing. 
It is written in C++ with Matlab_ interface. 
You can download it from github::

  git clone https://github.com/DiffusionMRITool/dmritool

Currently DMRITOOL_ has no GUI. It only provides executables for command line usage. 
You can also use the released mex executables in matlab_. 


With DMRITOOL_, you can 

* reconstruction/estimation of diffusion data, including diffusion weighted signal, ensemble average propagator (EAP), diffusion orientation distribution function (dODF), etc.
* generate spherically uniform sampling scheme design for single or multiple shells (will be in the future version).
* perform diffusion MRI data simulation.  
* visualization of spherical function field (e.g. dODF field, EAP profile field)



Citations will help us support the continued development of DMRITOOL_.
If you use the methods and codes released in DMRITOOL_, please cite the
related references. See `citation page <citation.html>`__.

Acknowledgements
================

DMRITOOL_ is/was supported by the following research groups:

*  `Section on Tissue Biophysics and Biomimetics (STBB)`_, `National Institutes of Health`_
*  `MIND Lab`_, `IDEA Group`_, `University of North Carolina at Chapel Hill`_
*  `Athena project team`_, INRIA_
*  `Brainnetome group`_, `Institute of Automation, Chinese Academy of Sciences`_


Users Documents
===============

.. toctree::
   :maxdepth: 2

   building
   userguide
   tutorials
   citation
   license

Development Documents
=====================

.. toctree::
   :maxdepth: 2

   developers
   codemodules


.. include:: links.inc
