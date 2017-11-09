% @MATLAB_FUNCTION_NAME@: get SPF basis matrix with given spherical and radial orders, orientationMatrix, bVector
%
% Usage:
%   spfMatrix = @MATLAB_FUNCTION_NAME@(shOrder, radialOrder, orientationMatrix, bVector)
%   spfMatrix = @MATLAB_FUNCTION_NAME@(shOrder, radialOrder, orientationMatrix, bVector, params)
% 
% Inputs:
%   shOrder:               SH order
%   radialOrder:           order in radial part
%   orientationMatrix :    orientation matrix where each row is an orientation (x,y,z), caresian format
%   bVector:               b values in a vector, each row is a b value
%   params.scale :         (optional) scale value for the SPF basis in radial part
%   params.tau:            (optional) tau value to calculate qVector from bVector, b=4\pi^2 \tau q^2
%   params.original:       (optional) true: (default) original SPF basis matrix 
%                                     false: independent SPF basis matrix, considering E(0)=1
%
% Outputs:
%   spfMatrix:              each column is a SPF basis
%
%
% Reference: 
%  1. "Model-Free, Regularized, Fast, and Robust Analytical Orientation Distribution Function Estimation",
%     Jian Cheng, Aurobrata Ghosh, Rachid Deriche, Tianzi Jiang, Medical Image Computing and Computer-Assisted Intervention (MICCAI'10), vol. 6361, pp. 648–656, sep, 2010.
%  2. "Model-free and Analytical EAP Reconstruction via Spherical Polar Fourier Diffusion MRI", 
%     Jian Cheng, Aurobrata Ghosh, Tianzi Jiang, Rachid Deriche, Medical Image Computing and Computer-Assisted Intervention (MICCAI'10), vol. 6361, pp. 590–597, sep, 2010. 
%
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%



