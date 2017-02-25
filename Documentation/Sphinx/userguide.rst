==========
User Guide
==========

.. meta::
   :description: user guide of dmritool
.. include:: meta_keywords.txt
.. include:: links.inc

.. contents:: Table of Contents
   :depth: 2
   :local:

4D image format and conversion
==============================

nifti format
------------
Diffusion weighted image (DWI) data is a 4D image data with 3D spatial dimensions and 1D diffusion dimension determined by the :math:`\mathbf{q}` value. 
In medical image data processing, we normally use nifti_ format to store the 4D image. 
.. See nifti header 

Multi-Volume image 
-------------------

The common way to store 4D image using nifti format is to sore different 3D volumes separately in memory. 
We call it as multi-volume image. 
FSL_ uses this way. 
When we use ITK_ to read this image format, we need to use `itk::Image\<double,4\> <http://www.itk.org/Doxygen/html/classitk_1_1Image.html>`__

In nifti format, a multi-volume image with size ``(64,64,20,120)`` is stored as ``dim[8] = 4 64 64 20 120 1 1 1``, where a scalar value is stored in each 4D voxel. 
See `dim description <http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/dim.html>`__.

VectorImage
-----------

In DMRITool_, we mainly use VectorImage format e.g. `ITK::VectorImage\<double,3\> <http://www.itk.org/Doxygen/html/classitk_1_1VectorImage.html>`__ for data processing. 
In this format, the values for the same voxel is stored continuously in memory. 
We can still nifti format to store it.  
We use this format because it is more efficient for processing 4D image data. 
See `this example <http://www.itk.org/Doxygen44/html/itkVectorImageTest_8cxx-example.html>`__. 

In nifti format, a vectorimage with size ``(64,64,20,120)`` is stored as ``dim[8] = 5 64 64 20 1 120 1 1``, where a 120 dimensional vector is stored in each 3D voxel. 
See `dim description <http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/dim.html>`__.


Conversion between multi-volume image and VectorImage
-----------------------------------------------------

We provide two executables for conversions between these two image storages. 

* ``VectorTo4DImageConverter`` is to convert a vectorimage to a multi-volume image.  ::

      VectorTo4DImageConverter  vectorimage.nii.gz  4dmultivolumeimage.nii.gz

* ``4DToVectorImageConverter`` is to convert a multi-volume image to a vectorimage. ::
      
      4DToVectorImageConverter  4dmultivolumeimage.nii.gz vectorimage.nii.gz  

* You can also use :doc:`ImageInfo <commands/ImageInfo>` to check the image header information ::
  
      ImageInfo vectorimage.nii.gz  
      ImageInfo 4dmultivolumeimage.nii.gz

* If you installed niftilib_, you can try::

      nifti_tool -disp_hdr -infiles vectorimage.nii.gz
      nifti_tool -disp_hdr -infiles 4dmultivolumeimage.nii.gz


DWI data input format
=====================

DWI data is a 4D image data where the fourth dimension is determined by a gradient file and a set of b values. 
It is possible to put gradient file and b values into the image header, 
as is done in NRRD_ format.  
Then the 3 files (e.g. image data file, gradient file, b values) becomes a single file. 

In DMRITool_, we choose to store the 3 files (e.g. image data file, gradient file, b values) separately, 
and the routines use a single text file which includes these 3 file names as input. 
The data reconstruction routines all rely on this input format. 

In this format, a 3 shell DWI data with b values of (1000,2000,3000) is stored with a single txt file (e.g. ``data.txt``). 
This txt file can have two different types. 

**Type 1**::

  ******in data.txt***********
  1000  grad1.txt  dwi1.nii.gz
  2000  grad2.txt  dwi2.nii.gz
  3000  grad3.txt  dwi3.nii.gz
  ****************************

Each line is for a single shell data. 
The gradient file ``grad1.txt`` stores a ``Nx3`` matrix. 
If the dimension ``N`` in ``grad1.txt`` is smaller than the dimension ``M`` in ``dwi1.nii.gz`` (dimension ``XxYxZxM``), 
the first ``M-N`` dimension in ``dwi.nii.gz`` is for b0 images.

.. Note::  DWI data files (e.g. ``dwi1.nii.gz``) can be a multi-volume image format or vectorimage format. Both formats are OK. 

If you just want to use some DWI volumes among these ``M`` DWI volumes, you can optionally add a index file in each line to specifically the index of volumes you want to use ::

  ******in data.txt**********************
  1000  grad1.txt  dwi1.nii.gz index1.txt
  2000  grad2.txt  dwi2.nii.gz index2.txt
  3000  grad3.txt  dwi3.nii.gz index3.txt
  ***************************************

The index file shows the index requested for reading. It starts from 0. 

**Type 2**: you can also specifically set all b values for all DWI volumes.  ::

 
  ******in data.txt*********************
  b1.txt  grad1.txt  dwi1.nii index1.txt
  b2.txt  grad2.txt  dwi2.nii index2.txt
  b3.txt  grad3.txt  dwi3.nii index3.txt
  **************************************

``b1.txt`` and ``grad1.txt`` should have the same dimension. The index files are optional. 
You can also catenate all b value files (gradients, DWI data files) into one b value file (gradients, DWI data files). ::

  ******in data.txt******************************
  ball.txt  gradall.txt  dwiall.nii indexall.txt
  ***********************************************

Configuration file format for DWI data simulation
===================================================

.. _DWIConfigurationFile:

We provides :doc:`DWISimulator <commands/DWISimulator>` routine to generate DWI data from a customizable configuration file.  
See example codes in ``Example/test.sh`` and example configuration files in ``Example`` folder. 

The configuration file specially sets the mixture of tensor model or mixture of configuration model for each voxel. 
The following table explains the parameters in the configuration file  

.. list-table::
   :header-rows: 1
   :widths: 30 30
   :stub-columns: 0
   
   *  -  DimSize
      -  output DWI image spatial size.

   *  -  ElementSpacing
      -  spatial spacing.

   *  -  Scale
      -  DWI image for ``b=0``.

   *  -  ModelType
      -  Symmetric tensor model, or cylinder model, spherical coordinate or Cartesian coordinate.

   *  -  DiffusionParameters
      -  First three numbers are for spatial position. The following numbers are respectively for each tensor/cylinder components. |br|
         For example, when ``ModelType=SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS``, |br|
         ``DiffusionParameters = 3 3 0 90 0 0.0015 0.0003 0.5 90 60 0.0015 0.0003 0.5`` means in voxel ``(3,3,0)``, |br| 
         two tensors with eigenvalues ``(0.0015,0.0003)``, partial volume weight ``0.5``, spherical coordinates (90,0) and (90,60) respectively. |br|
         Thus crossing angle is ``60`` degree. 

   *  -  BackgroundDiffusionParameters
      -  diffusion parameters for all other voxels which are not specially set in ``DiffusionParameters``.

   *  -  RicianNoiseSigma
      -  Sigma for  Rician noise. Suggest to set it as 0, then use :doc:`DWINoiseGenerator <commands/DWINoiseGenerator>` to add noise to the noise-free DWI data.




Visualization of spherical function field
==========================================

Generate and visualize vtk files 
--------------------------------

Both ODFs and EAP profile are spherical function fields, where in each voxel there is a spherical function (e.g. ODF, EAP profile with a fixed radius). 
The spherical function field can be represented using uniform samples of spherical functions or the Spherical Harmonics (SH) coefficients of spherical functions. 

We provides some routines to generate mesh files (.vtk format) from the SH coefficients or spherical function samples. 
Then these vtk files can be visualized in paraview_ 
or using the routine :doc:`vtkviewer <commands/vtkviewer>` 
and :doc:`VTKPolyData.py <commands/VTKPolyData.py>`. 
Please see :ref:`the command list for visualization <commandlist_Visualization>`. 
To use :doc:`VTKPolyData.py <commands/VTKPolyData.py>`, you only need to build VTK with python wrapping. 
It does not need to build dmritool_. 

Visualization using Paraview
----------------------------

When using paraview_, it is possible to visualize these generated vtk files with a scalar map, e.g. the GFA map, as the background. 
Paraview_ originally can visualize vtk files.
For the GFA map which is in nifti format, paraview can also visualize it when the ``AnalyzeNiftiIO`` plugin is enabled.
You can enable it by clicking ``Tools`` -> ``Manage Plugins`` -> ``AnalyzeNiftiIO``.  

Note that it supports ``.nii`` file, but not ``.nii.gz`` file, and you may need to set the origin manually in paraview such that
the gfa map can be visualized in the same coordinate as the EAP profile in vtk format.

* In the ``properties`` panel, when you click ``toggle advanced properties``, it shows ``Translation`` and ``Scale``.
* By setting ``translation``, you can translate the image map to different position.
* By setting ``scale``, you scale the axis. If you set ``(-1,1,1)`` for scale, it reverses the x-axis.

You can also enable other useful plugins for paraview, e.g. ``quadview``.
 



