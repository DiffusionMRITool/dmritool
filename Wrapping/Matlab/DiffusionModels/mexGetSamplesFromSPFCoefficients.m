% @MATLAB_FUNCTION_NAME@: get DWI or EAP samples (orignal samples or the SH coefficients) from given SPF coefficients
%
% Usage:
%   [samples] = @MATLAB_FUNCTION_NAME@(spfCoef4D, radius, params)
%   [samples] = @MATLAB_FUNCTION_NAME@(spfCoef4D, radiusVec, grad, params)
% 
% Inputs:
%   dwi4DNoramlized:       4D DWI data, normalized by b0 image;
%   grad :                 orientation matrix where each row is an orientation (x,y,z), caresian format. 
%                          If it is not given, then return SH coefficients instead of orignal samples.
%   radiusVector:          b values in a vector, or r values. If it is set, the number of rows should be the same as grad
%   radius:                b value or r value. 
%
%   params.MD0 :           Typical mean diffusivity (MD) value. The default value is 0.7e-3. 
%   params.tau :           Tau value. The default is calculated based on 4*pi*pi*tau=1.
%   param.scale :          Scale value for the SPF basis in radial part, the default value is set by tau and MD0 and tau 
%   params.sh  :           SH order
%   params.ra  :           order in radial part
%   params.fourier :       Use Fourier transform. 
%   params.inqspace :      If it is true, the basis is in q-space (default). If it is false, the basis is in r-space. 
%   params.mdImage  :      Mean diffusivity Image for adaptive scale. 
%   params.basisType :     SPF: use SPF basis (default). 
%   params.mask :          Mask file. 2D or 3D matrix. 
%
%   params.thread :        Number of thread. Default is -1, which means it is automatically determined.
%   params.verbose :       Verbose level. 0: no output log. 1: normal log. 2: large log for debug.
%
% Outputs:
%   isoDWI:                each column is a dwi sample
%
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%

