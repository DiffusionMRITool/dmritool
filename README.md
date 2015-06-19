DMRITOOL 
========

[![Build Status](https://travis-ci.org/DiffusionMRITool/dmritool.svg?branch=master)](https://travis-ci.org/DiffusionMRITool/dmritool)

Introduction
============

[DMRITOOL](http://diffusionmritool.github.io/) is a **free** and **open source** toolbox for [diffusion MRI](http://en.wikipedia.org/wiki/Diffusion_MRI) data processing. 
It is written in C++ with matlab interface. 

With DMRITOOL, you can:

* perform reconstruction/estimation of diffusion data, including diffusion weighted signal, ensemble average propagator (EAP), diffusion orientation distribution function (dODF), and some meaningful scalar maps, etc.
* generate spherically uniform sampling schemes for single or multiple shells.
* perform diffusion MRI data simulation.
* visualize spherical function field (e.g. dODF field, EAP profile field)



Website
=======

Please check the [DMRITOOL website](http://diffusionmritool.github.io/) for documentation and more information.

Download
========

You can download the latest source codes from [github](https://github.com/DiffusionMRITool/dmritool):

    git clone https://github.com/DiffusionMRITool/dmritool

Building
========

See [this page](http://diffusionmritool.github.io/building.html) for building the source codes. 

Citation
========

Citations will help us support the continued development of DMRITOOL. 

If you use the methods and codes released in DMRITOOL, please cite the related references. 
See [the citation page](http://diffusionmritool.github.io/citation.html) for details. 

Acknowledgements
================

DMRITOOL is/was supported by the following research groups:

* [Section on Tissue Biophysics and Biomimetics (STBB)](http://stbb.nichd.nih.gov/index.html), [National Institutes of Health](http://www.nih.gov/)
* [MIND Lab](http://www.unc.edu/~ptyap/index.html), [IDEA Group](https://www.med.unc.edu/bric/ideagroup), [University of North Carolina at Chapel Hill](http://www.unc.edu/)
* [Athena project team](https://team.inria.fr/athena/), [INRIA](http://www.inria.fr/)
* [Brainnetome center](http://www.brainnetome.org/en/), [Institute of Automation, Chinese Academy of Sciences](http://english.ia.cas.cn/)

License
=======

DMRITOOL is a free open source software. 
It is currently under the [GNU General Public License](http://www.gnu.org/licenses/gpl.html), 
because it uses [GSL](http://www.gnu.org/software/gsl/) for mathematical special functions and [SPAMS](http://spams-devel.gforge.inria.fr/) for some optimization problems. 

The software under the license is distributed on an "as is" basis, without warranties.
It is user's responsibility to validate the behavior of the routines and their accuracy using the released source code. 

