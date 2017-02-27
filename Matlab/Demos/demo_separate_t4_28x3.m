%% Uniformly select multiple subsets from a given set
%
% This is a demo to uniformly separate multiple subsets from a given set by 
% using Mixed Integer Linear Programming (MILP)
% It reproduces the experiment in the paper. 
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

%% extract 28x3 samples from grad_t4 using MILP
clear params grbParams
params.numSamples = [28, 28, 28];
% sos constraint
params.sos = 1;
% set a lower bound, it can be 0 or covering radius from an existing scheme
params.lbCost = [0.38,0.38,0.38,0.23];


params.w=0.5;

% grb parameters
% MIPFocus=1 seems better
grbParams.MIPFocus = 1;
% time limit
grbParams.TimeLimit = 600;
% print verbose output from gurobi
grbParams.OutputFlag=true;
% grbParams.OutputFlag=false;

% It is possible to set a warm start if you have one
% grbParams.start = gradInit;

% run
[grad,grb, indexMatrix] = OptimalSamplingMultiSubsetsFromSameSet(grad_t4, params,grbParams);

%% visualize the result
VisualizeMultiShellScheme(grad{1},grad{2},grad{3});
title({'Estimated Scheme.', ['Combined covering radius = ', num2str(CoveringRadius([grad{1};grad{2};grad{3}])*180/pi), ' degree']});
