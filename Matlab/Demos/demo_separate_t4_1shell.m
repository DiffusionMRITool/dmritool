%% Uniformly select a single subset from a given set
%
% This is a demo to uniformly separate a subset from a given set by 
% using Mixed Integer Linear Programming (MILP)
% 
% Reference:
% "Designing Single- and Multiple-Shell Sampling Schemes for Diffusion MRI Using Spherical Code", 
% Jian Cheng, Dinggang Shen, Pew-Thian Yap, MICCAI 2014.  
%
% Copyright (c) 2013, Jian Cheng (jian.cheng.1983@gmail.com)
%

%% read sphere tessellation with 321 samples in the hemisphere
grad_t4 = ReadDirections([getenv('HOME'), '/.dmritool/Data/Tessellation/directions_t4.txt']);

VisualizeMultiShellScheme(grad_t4);
title(['Sphere tessellation with 321 samples in the hemisphere']);

%% extract 30 samples from grad_t4 using MILP
clear params grbParams
params.numSamples = 30;
% set a lower bound, it can be 0 or covering radius from an existing scheme
params.lbCost = 0.24;

% grb parameters
% MIPFocus 1 seems better
grbParams.MIPFocus=1;
% set time limit as 10 minutes or more
grbParams.TimeLimit=600;
% print verbose output from gurobi
grbParams.OutputFlag=true;
% grbParams.OutputFlag=false;
% params.ModelFile='model.mps';

% run
[grad,grb] = OptimalSamplingSingleSubset(grad_t4, params, grbParams);

%% visualize the result
VisualizeMultiShellScheme(grad);
title(['Estimated Scheme. Covering radius = ', num2str(CoveringRadius(grad)*180/pi), ' degree'], 'FontSize', 15);
