%% Generate a multi-shell scheme by using IMOC + 1-Opt + CNLO
%
% This is a demo to generate a multi-shell scheme (e.g. 3 shells 28 samples per shell)
% by using IMOC + 1-Opt + CNLO. 
% 
% OptimalSamplingMultiShellCNLO is in matlab. But to run IMOC and 1-Opt, you need to build the dmritool C++ codes.
% 
% Reference:
%   "Novel single and multiple shell uniform sampling schemes for diffusion MRI using spherical codes", 
%   Jian Cheng, Dinggang Shen, Pew-Thian Yap, Peter J. Basser, MICCAI 2015.  
%
% Copyright (c) 2016, Jian Cheng (jian.cheng.1983@gmail.com)
%

%% IMOC to design a multi-shell scheme (28,28,28)
tic
! SamplingSchemeQSpaceIMOCEstimation grad_28x3_IMOC.txt --tessOrder 7 --numberOfSamples 28,28,28 
toc

grad_IMOC_shell1 = ReadDirections('grad_28x3_IMOC_shell1.txt');
grad_IMOC_shell2 = ReadDirections('grad_28x3_IMOC_shell2.txt');
grad_IMOC_shell3 = ReadDirections('grad_28x3_IMOC_shell3.txt');

fprintf('\nConvering radii:\t %10.3f\t%10.3f\t%10.3f\t%10.3f\n',...
    CoveringRadius(grad_IMOC_shell1)*180/pi,CoveringRadius(grad_IMOC_shell2)*180/pi,CoveringRadius(grad_IMOC_shell3)*180/pi,...
    CoveringRadius([grad_IMOC_shell1;grad_IMOC_shell2;grad_IMOC_shell3])*180/pi);

VisualizeMultiShellScheme(grad_IMOC_shell1, grad_IMOC_shell2, grad_IMOC_shell3);
title({'Scheme by IMOC .', ['Combined covering radius = ', num2str(CoveringRadius([grad_IMOC_shell1;grad_IMOC_shell2;grad_IMOC_shell3])*180/pi), ' degree']});



%% IMOC + 1-Opt to design a multi-shell scheme (28,28,28)
tic
! SamplingSchemeQSpace1OptEstimation grad_28x3_IMOC1Opt.txt --initial grad_28x3_IMOC_shell1.txt,grad_28x3_IMOC_shell2.txt,grad_28x3_IMOC_shell3.txt --tessOrder 7 
toc

grad_IMOC_1Opt_shell1 = ReadDirections('grad_28x3_IMOC1Opt_shell1.txt');
grad_IMOC_1Opt_shell2 = ReadDirections('grad_28x3_IMOC1Opt_shell2.txt');
grad_IMOC_1Opt_shell3 = ReadDirections('grad_28x3_IMOC1Opt_shell3.txt');

fprintf('\nConvering radii:\t %10.3f\t%10.3f\t%10.3f\t%10.3f\n',...
    CoveringRadius(grad_IMOC_1Opt_shell1)*180/pi,CoveringRadius(grad_IMOC_1Opt_shell2)*180/pi,CoveringRadius(grad_IMOC_1Opt_shell3)*180/pi, ...
    CoveringRadius([grad_IMOC_1Opt_shell1;grad_IMOC_1Opt_shell2;grad_IMOC_1Opt_shell3])*180/pi);

VisualizeMultiShellScheme(grad_IMOC_1Opt_shell1, grad_IMOC_1Opt_shell2, grad_IMOC_1Opt_shell3);
title({'Scheme by IMOC + 1-Opt.', ['Combined covering radius = ', num2str(CoveringRadius([grad_IMOC_1Opt_shell1;grad_IMOC_1Opt_shell2;grad_IMOC_1Opt_shell3])*180/pi), ' degree']});


%% IMOC + 1-Opt + CNLO to design a multi-shell scheme (28,28,28)

clear param;

% weight between individual shell and combined shell
param.w  = 0.5;

% solver can be nlopt or matlab. Matlab solver is slow.
param.solver  = 'nlopt';
% param.solver  = 'matlab';

% use local constraint for fast CNLO
param.localCon = 2;
param.localConAngle = 0.1;

% maxtime only for single run, the run time may be longer for whole process if the number of samples is large.
param.maxtime = 600;
param.verbose = 1;

grad{1}=grad_IMOC_1Opt_shell1;
grad{2}=grad_IMOC_1Opt_shell2;
grad{3}=grad_IMOC_1Opt_shell3;

tic
[gradCell,xopt, fopt, retcode] = OptimalSamplingMultiShellCNLO(grad, param);
toc

fprintf('\nConvering radii:\t %10.3f\t%10.3f\t%10.3f\t%10.3f\t\n',...
    CoveringRadius(gradCell{1})*180/pi,CoveringRadius(gradCell{2})*180/pi,CoveringRadius(gradCell{3})*180/pi, ...
    CoveringRadius([gradCell{1};gradCell{2};gradCell{3}])*180/pi);

VisualizeMultiShellScheme(gradCell{1}, gradCell{2}, gradCell{3});
title({'Scheme by IMOC + 1-Opt + CNLO.', ['Combined covering radius = ', num2str(CoveringRadius([gradCell{1};gradCell{2};gradCell{3}])*180/pi), ' degree']});

