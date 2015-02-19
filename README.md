DMRITOOL 
=========

[![Build Status](https://travis-ci.org/DiffusionMRITool/dmritool.svg?branch=master)](https://travis-ci.org/DiffusionMRITool/dmritool)
[![Documentation Status](https://readthedocs.org/projects/dmritool/badge/?version=latest)](https://readthedocs.org/projects/dmritool/?badge=latest)

Introduction
============

[DMRITOOL](http://dmritool.readthedocs.org/en/latest/index.html) is a **free** and **open source** toolbox for [diffusion MRI](http://en.wikipedia.org/wiki/Diffusion_MRI) data processing. 
It is written in C++ with matlab interface. 

DMRITOOL mainly focus on:

* Reconstruction/estimation of diffusion data, 
  including diffusion weighted signal, ensemble average propagator (EAP), diffusion orientation distribution function (dODF), 
  and some meaningful scalar maps. 

* Spherically uniform sampling scheme design for single or multiple shells. 

Documentation
=============

Please check the [documentation](http://dmritool.readthedocs.org/en/latest/index.html).


Building
========

See [this page](http://dmritool.readthedocs.org/en/latest/building.html) for building the source codes. 

Citation
========

Citations will help us support the continued development of DMRITOOL. 

If you use the methods and codes released in DMRITOOL, please cite the related references. 
See [the citation page](http://dmritool.readthedocs.org/en/latest/citation.html) for details. 

Acknowledgements
================

DMRITOOL is/was supported by the following research groups:

* [Section on Tissue Biophysics and Biomimetics (STBB)](http://stbb.nichd.nih.gov/index.html), [National Institutes of Health](http://www.nih.gov/)
* [MIND Lab](http://www.unc.edu/~ptyap/index.html), [IDEA Group](https://www.med.unc.edu/bric/ideagroup), [University of North Carolina at Chapel Hill](http://www.unc.edu/)
* [Athena project team](https://team.inria.fr/athena/), [INRIA](https://team.inria.fr/athena/)
* [Brainnetome group](http://www.brainnetome.org/en/), [Institute of Automation, Chinese Academy of Sciences](http://english.ia.cas.cn/)

License
=======

DMRITOOL is a free open source software. 
It is currently under the [GNU General Public License](http://www.gnu.org/licenses/gpl.html), 
because it uses [GSL](http://www.gnu.org/software/gsl/) for mathematical special functions and [SPAMS](http://spams-devel.gforge.inria.fr/) for some optimization problems. 

The software under the license is distributed on an "as is" basis, without warranties.
It is user's responsibility to validate the behavior of the routines and their accuracy using the released source code. 

