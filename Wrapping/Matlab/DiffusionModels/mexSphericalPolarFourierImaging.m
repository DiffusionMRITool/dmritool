% @MATLAB_FUNCTION_NAME@: estimate continuous DWI signal and EAP represented by coefficients of given basis. 
%
% Usage:
%   [signal] = @MATLAB_FUNCTION_NAME@(dwi4DNoramlized, grad, bVec, params)
% 
% Inputs:
%   dwi4DNoramlized:       4D DWI data, normalized by b0 image;
%   grad :                 orientation matrix where each row is an orientation (x,y,z), caresian format
%   bVector:               b values in a vector, each row is a b value
%
%   params.MD0 :           Typical mean diffusivity (MD) value. The default value is 0.7e-3. 
%   params.tau :           Tau value. The default is calculated based on 4*pi*pi*tau=1.
%   param.scale :          Scale value for the SPF basis in radial part, the default value is set by tau and MD0 and tau 
%   params.sh:             SH order
%   params.ra:             order in radial part
%   params.lambdaSH :      lambda in spherical part
%   params.lambdaRA :      lambda in radial part
%   params.lambdaL1 :      lambda for L1 solver
%   params.odfOrder :      odf order 
%   params.radius   :      Radius for EAP profile.
%   params.mdImage  :      Mean diffusivity Image for adaptive scale.
%   params.estimation :    (LS, L1_2, L1_DL) LS: least square. 
%                          L1_2 means Laplace-Beltrami weight (two regularization parameters lambdaSH and lambdaRA); 
%                          L1_DL only use lambdaL1, because the learned dictionary corresponses the regularization matrix.
%   params.solver   :      (FISTA_LS, SPAMS)
%                          FISTA_LS: FISTA using least square initialization.  
%                          SPAMS (default): use spams' weighted lasso solver. 
%   params.maxIter  :      Maximal number of iteration in L1 FISTA.
%   params.minChange :     Minimal change percentage of the cost function and variable for l1 oprimization.
%   params.mask :          Mask file. 2D or 3D matrix. 
%
%   params.thread :        Number of thread. Default is -1, which means it is automatically determined.
%   params.verbose :       Verbose level. 0: no output log. 1: normal log. 2: large log for debug.
%
% Outputs:
%   signal:                coefficient under a given basis, 4D image
%
%
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%


