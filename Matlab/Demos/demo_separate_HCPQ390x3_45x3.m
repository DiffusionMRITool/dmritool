%% Uniformly select multi-shell subsets from a multi-shell scheme
%
% This is a demo to uniformly separate multiple subsets from the multi-shell scheme in
% Human Connectome Project (HCP) using Mixed Integer Linear Programming (MILP)
% 
%
% Reference:
%
% 1. "Single- and Multiple-Shell Uniform Sampling Schemes for Diffusion MRI Using Spherical Codes", 
%     Jian Cheng, Dinggang Shen, Pew-Thian Yap, Peter J. Basser, IEEE Transactions on Medical Imaging, 2017.  
%
% 2. "Designing Single- and Multiple-Shell Sampling Schemes for Diffusion MRI Using Spherical Code", 
%     Jian Cheng, Dinggang Shen, Pew-Thian Yap, MICCAI 2014.  
%
% Copyright (c) 2013, Jian Cheng (jian.cheng.1983@gmail.com)
%

%% read HCP scheme
grad_hcp_all{1}=ReadDirections('/home/jcheng/mywork/output/sampling/HCP_Q3/grad_b1000.txt');
grad_hcp_all{2}=ReadDirections('/home/jcheng/mywork/output/sampling/HCP_Q3/grad_b2000.txt');
grad_hcp_all{3}=ReadDirections('/home/jcheng/mywork/output/sampling/HCP_Q3/grad_b3000.txt');

VisualizeMultiShellScheme(grad_hcp_all{1}, grad_hcp_all{2}, grad_hcp_all{3});
title(['Original 90x3 HCP scheme']);

%% extract 45x3 samples from 90x3 scheme
clear params grbParams
params.numSamples = [45, 45, 45];
% set a lower bound, it can be 0 or covering radius from an initial guess
params.lbCost = [0.25,0.25,0.25,0.11];
% params.lbCost = [0,0,0,0];

% if you want to write the model for gurobi_cl 
% params.ModelFile='model.mps'
% params.ResultFile='grb_solution.sol'


params.w=0.5;

% grb parameters
% MIPFocus 1
grbParams.MIPFocus=1;
grbParams.TimeLimit=600;
grbParams.OutputFlag=true;


% run
[grad,grb, indexMatrix] = OptimalSamplingMultiSubsetsFromDifferentSets(grad_hcp_all, params, grbParams);

% the selected indices for each shell
index_1=find(indexMatrix{1}==1);
index_2=find(indexMatrix{2}==1);
index_3=find(indexMatrix{3}==1);

%% visualize the result
VisualizeMultiShellScheme(grad{1},grad{2},grad{3});
title({'Estimated scheme.', ['Combined covering radius = ', num2str(CoveringRadius([grad{1};grad{2};grad{3}])*180/pi), ' degree']});

