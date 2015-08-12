%% Uniformly separate several subsets from a given set
%
% This is a demo to uniformly separate two subsets from a given set by 
% using Mixed Integer Linear Programming (MILP). 
% It reproduces the experiment in the paper. 
%
% Reference:
% "Designing Single- and Multiple-Shell Sampling Schemes for Diffusion MRI Using Spherical Code", 
% Jian Cheng, Dinggang Shen, Pew-Thian Yap, MICCAI 2014.  
%
% Copyright (c) 2013, Jian Cheng (jian.cheng.1983@gmail.com)
%

%% Read two sets of uniform directions
grad1 = ReadDirections([getenv('HOME'), '/.dmritool/Data/Tessellation/directions_t3.txt']);
n1 = size(grad1, 1);

grad2 = ReadDirections([getenv('HOME'), '/.dmritool/Data/ElectricRepulsion/Elec060.txt']);
n2 = size(grad2, 1);

% Visualize the original two schemes. Different colors denote samples in different shells.
VisualizeMultiShellScheme(grad1, grad2);
title(['Orignal two schemes']);

%% Randomly mix these two sets
index = randperm(n1+n2);

gradAll = [grad1; grad2];
gradAll = gradAll(index,:);

% Visualize the mixed two schemes. 
VisualizeMultiShellScheme(gradAll);
title(['Randomly mixed scheme']);

%% Extract two sets using MILP
% Now we want to seperate these two sets from the mixture. 

clear params grbParams
% set parameters
% w=1 means we do not care about the combined shell with all samples. See the paper. 
param.w=1;
param.numSamples=[n1 n2];

% use default grb parameters
[gradCell, grb, indexMatrix]=OptimalSamplingMultiSubsetsFromSameSet(gradAll,param, []);

%% Test the results 
% covering radius comparison
fprintf('covering radius of set 1 = %f degree\n', CoveringRadius(grad1)*180/pi);
fprintf('covering radius of set 2 = %f degree\n', CoveringRadius(grad2)*180/pi);
fprintf('covering radius of extracted set 1 = %f degree\n', CoveringRadius(gradCell{1})*180/pi);
fprintf('covering radius of extracted set 2 = %f degree\n', CoveringRadius(gradCell{2})*180/pi);

% number of wrongly detected samples
fprintf('the number of extracted directions not in set 1: %d\n', n1-sum(indexMatrix(index<=n1,1)))
fprintf('the number of extracted directions not in set 2: %d\n', n2-sum(indexMatrix(index>n1,2)))

% Visualize the separated two sets
VisualizeMultiShellScheme(gradCell{1}, gradCell{2});
title(['Separated two schemes']);

