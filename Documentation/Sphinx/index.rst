========================
DMRITool's documentation
========================

.. meta::
   :description: homepage of dmritool
.. include:: meta_keywords.txt
.. include:: links.inc

.. contents:: Table of Contents
   :depth: 2
   :local:

.. image:: https://travis-ci.org/DiffusionMRITool/dmritool.svg?branch=master
    :target: https://travis-ci.org/DiffusionMRITool/dmritool
    :alt: travis-ci Status             

.. image:: https://api.codacy.com/project/badge/Grade/13b3fe25b4b84ac6a413bb31e29d58fc    
    :target: https://www.codacy.com/app/JianCheng/dmritool?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DiffusionMRITool/dmritool&amp;utm_campaign=Badge_Grade
    :alt: Codacy Status             


Introduction
============

DMRITool_ is a **free** and **open source** toolbox for `diffusion MRI <http://en.wikipedia.org/wiki/Diffusion_MRI>`__ data processing. 
It is written in C++ with Matlab_ interface. 
Currently DMRITool_ has no GUI. It only provides executables for command line usage. 
You can also use the released mex executables in matlab_. 


With DMRITool_, you can: 

* perform reconstruction/estimation of diffusion data, including diffusion weighted signal, 
  ensemble average propagator (EAP), diffusion orientation distribution function (dODF), 
  and some meaningful scalar maps, etc.
* generate spherically uniform sampling schemes for single or multiple shells.
* perform diffusion MRI data simulation.  
* visualize spherical function field (e.g. dODF field, EAP profile field)


Citations will help us support the continued development of DMRITool_.
If you use the methods and codes released in DMRITool_, please cite the
related references. See :doc:`citation page <citation>`.

Download
========

You can download the latest source codes from `github <https://github.com/DiffusionMRITool/dmritool>`__:


.. raw:: html
 
  <div>
    <IMG SRC="_static/GitHub.png" height="40" width="40"> 
    <a href="https://github.com/DiffusionMRITool/dmritool/zipball/master"> [ZIP]   </a>  
    <a href="https://github.com/DiffusionMRITool/dmritool/tarball/master">[TAR.GZ] </a>
    <p></p>
  </div>

Or use git_:

:: 

  git clone https://github.com/DiffusionMRITool/dmritool



Acknowledgements
================

DMRITool_ is/was supported by the following research groups:

*  `Section on Quantitative Imaging and Tissue Sciences (SQITS)`_, `National Institutes of Health`_
*  `MIND Lab`_, `IDEA Group`_, `University of North Carolina at Chapel Hill`_
*  `Athena project team`_, INRIA_
*  `Brainnetome center`_, `Institute of Automation, Chinese Academy of Sciences`_


Users Documentation
===================

.. sidebar:: Contact Us 
   
   If you have questions about dmritool, please contact us: dmritool-discussion@www.nitrc.org

.. toctree::
   :maxdepth: 2

   Home <self>
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
   Doxygen Documentation <http://diffusionmritool.github.io/dmritool-doxygen>
   codemodules


