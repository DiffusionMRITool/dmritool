function [gradCell, xopt, fopt, retcode] = OptimalSamplingMultiShellCNLO(gradCellInitial, param)
% update gradients from an initial gradient set, such that the updated gradients are evenly distributed. 
%
% USAGE:
%   [gradCell, xopt, fopt, retcode] = OptimalSamplingMultiShellCNLO(gradCellInitial, param)
%
% INPUT (required)
%   gradCellInitial      :  N_k x 3 gradient matrices, where each row is a point in sphere. k=1,2,...,K
%   param.w              :  weight for single shell term. 0<w<1. 0.5 is the default value.
% 
% INPUT (optional)
%   param.solver         :  'nlopt' (default): use sqp solver (SLSQP) in NLOPT.
%                           'matlab': use fmincon (sqp solver).
%   param.cartesian      :  true (default): use cartesian coordinate (with unit equality constraint);
%                           false: use spherical coordinate (no unit equality constraint).
% 
%   param.localCon       :  0: do not use local constraint
%                           1: use local constraint for the distance between samples in results and samples in initialization to reduce the inequality constraint. 
%                              It is faster, and with an appropriate localConAngle, it can obtain the same results as it is set as false. 
%                           2 (default): repeat the CNLO with a fixed localConAngle until the result keeps the same.
%
%   param.localConAngle  :  It is the maximal angular change between samples in results and samples in initialization. 
%                           It is used when param.localCon is 1 or 2. The default value is 0.1 radian, if it is not given. 
%                           If localConAngle>pi/2, then the result is the same as the results without localCon.
%
% 
% OUTPUT
%   gradCell             :  Kx1 cell, each element is a gradient matrix. 
%
%
% Reference:
%   1. "Single- and Multiple-Shell Uniform Sampling Schemes for Diffusion MRI Using Spherical Codes", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, Peter J. Basser, IEEE Transactions on Medical Imaging, 2017.  
%   2. "Novel single and multiple shell uniform sampling schemes for diffusion MRI using spherical codes", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, Peter J. Basser, MICCAI 2015.  
%
% Copyright (c) 2015, Jian Cheng <jian.cheng.1983@gmail.com>

%% sampling parameters

paramDefault.w=0.5;
paramDefault.solver='nlopt';
paramDefault.cartesian=true;
paramDefault.localCon=2;
paramDefault.localConAngle=0.1;

paramDefault.maxtime=-1;
paramDefault.verbose=0;

if nargin==1
    param = paramDefault;
else
    param = CopyDefaultStructParams(paramDefault, param);
end

paramInner = param;

if param.localCon==2
    paramInner.localCon = 1;
end

%%
numShells = numel(gradCellInitial);
numSamples = zeros(numShells,1);
for i = 1 : numShells
    [numSamples(i),c] = size(gradCellInitial{i});
    if c~=3
        error('wrong size of input gradients!');
    end
end

upperBoundSpherical = CoveringRadiusUpperBound(2*numSamples);
upperBoundSphericalAll = CoveringRadiusUpperBound(2*sum(numSamples));

if param.localCon>=1
    if ~isfield(param, 'localConAngle')
        if param.localCon==2
%             param.localConAngle = min([upperBoundSphericalAll; upperBoundSpherical]) / 2.0;
            param.localConAngle = 0.1;
        else
            param.localConAngle = max([upperBoundSphericalAll; upperBoundSpherical]) / 1.5;
        end
    end
end
paramInner.localConAngle = param.localConAngle;

%% CNLO loop


if (param.localCon==2)
    
    isStable = false;
    gradCellOld = gradCellInitial;
    gradAllOld = GetGradAllShells(gradCellOld);

    while(~isStable)
        
        [gradCell, xopt, fopt, retcode] = OptimalSamplingMultiShellCNLO_singleRun(gradCellOld, paramInner);
        
        isStable = true;
        for i = 1:numShells
            angle = CoveringRadius(gradCell{i})*180/pi;
            angleOld = CoveringRadius(gradCellOld{i})*180/pi;;
            if (param.verbose)
                fprintf('%d-th shell: \tangleOld=%10.10f\tangle=%10.10f\n', i, angleOld, angle)
            end
            if (angle-angleOld)/angleOld > 1e-5
                isStable = false;
                break;
            end
        end
        gradAll = GetGradAllShells(gradCell);
        angle = CoveringRadius(gradAll)*180/pi;
        angleOld = CoveringRadius(gradAllOld)*180/pi;
        if (param.verbose)
            fprintf('Combined shell: \tangleOld=%10.10f\tangle=%10.10f\n', angleOld, angle)
        end
        if (angle-angleOld)/angleOld > 1e-5
            isStable = false;
        end
        
        gradCellOld = gradCell;
        gradAllOld = gradAll;
    end
else 
     [gradCell, xopt, fopt, retcode] = OptimalSamplingMultiShellCNLO_singleRun(gradCellInitial, paramInner);
end


%% output
index_shells=0;
for k = 1 : numShells
    if param.cartesian
        grad_i = xopt(numShells+1+index_shells+1 : numShells+1+index_shells + 3*numSamples(k));
        gradCell{k} = reshape(grad_i, [3,numSamples(k)])';
        index_shells = index_shells + 3*numSamples(k);
    else
        grad_i = xopt(numShells+1+index_shells+1 : numShells+1+index_shells + 2*numSamples(k));
        grad_i = reshape(grad_i, [2,numSamples(k)])';
        [gradCell{k}(:,1), gradCell{k}(:,2), gradCell{k}(:,3)] = SphericalToCartesian(grad_i(:,1), grad_i(:,2), 1);
        index_shells = index_shells + 2*numSamples(k);
    end
end

