% @MATLAB_FUNCTION_NAME@: write Matlab image matrix to any image format supported by ITK using itk::Image and itk::VectorImage.
%
% Usage:
%   @MATLAB_FUNCTION_NAME@('filename', image, [params])
%
% INPUT
%  filename   :     file name to save. 
%  image      :     4D matlab matrix
%
%  params.origin            :  origin in image 
%  params.spacing           :  spacing in image 
%  params.vectorImage       :  If it is true, save 3D itk::VectorImage, otherwise save 4D itk::Image (default). 
% 
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%
