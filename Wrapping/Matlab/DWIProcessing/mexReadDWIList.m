% @MATLAB_FUNCTION_NAME@: read dwi images, gradient files, and b values
%
% Usage:
%   [dwi4D, bVector, gradMatrix, b0Image3D] = @MATLAB_FUNCTION_NAME@(txtFile, params)
%   [dwi4D, bVector, gradMatrix] = @MATLAB_FUNCTION_NAME@(txtFile)
%
% INPUT
%   txtFile    : file name which includes DWI files, b values, gradient files. 
%                * ******in data.txt*******************
%                * 650   grad1.txt  dwi1.nii index1.txt
%                * 1500  grad2.txt  dwi2.nii index2.txt
%                * 3000  grad3.txt  dwi3.nii index3.txt
%                * ************************************
%                * Each line is for a single shell data. 
%                * If the dimension N in grad.txt is smaller than the dimension M in dwi.nii, the first M-N dimension in dwi.hpp is for b0 image
%                * index.txt (optional) shows the index requested for reading.  
%                * dwi.hdr can be stored as VectorImage or Image 
%                *
%                * ******in data.txt*********************
%                * b1.txt  grad1.txt  dwi1.nii index1.txt
%                * b2.txt  grad2.txt  dwi2.nii index2.txt
%                * b3.txt  grad3.txt  dwi3.nii index3.txt
%                * **************************************
%                * b.txt and grad.txt should have the same dimension. In this format, it is not necessary to have the same b values in b1.txt, b2.txt, etc
%
%
%   params.b0Image    :    A 3D image for b0. If it is set, then the b0 image should not be in the data.txt, 
%                          which means the dimension of grad1.txt and dwi.nii.gz should have the same size. 
%                          If it is not set and grad1.txt and dwi.nii.gz have the same size, then the values of b0Image are all one.
%
%   params.normalize  :    If it is true, then normalize the DWI data using b0 image. 
%                          The default value is true
%
%   params.bThreshold :    If it is positive, then correct b values by grouping the b values whose difference is smaller than the bThreshold. 
%                          The default value is -1
%
%   params.correctDWI :    correct DWI values if it is more than the value in b0 image or if it is no more than 0.
%   params.warn :          if it is true, print some warnings if the dwi values are abnormal (i.e. zero, negative, larger than b0 value)
%
%              
%
% OUTPUT
%   dwi4D   :  4D dwi image
%   bVector :  vector of b values
%   grad    :  Nx3 matrix, each row is a point in sphere.
%   b0Image :  3D b0 image
%
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%

