function [gradCell,grb_result,indexMatrix] = OptimalSamplingMultiSubsetsFromDifferentSets(gradAll, param, grbParam)
% get K subsets of gradients from corresponding K given gradient sets, such that
% gradients in each subset are evenly distributed. 
%
% USAGE:
%   [gradCell,grb_result,indexMatrix] = OptimalSamplingMultiSubsetsFromDifferentSets(gradAll, param, grbParam)
%
% INPUT
%   gradAll            :  K gradient matrix cell, where each element is a N_kx3 gradient matrix for a single shell.
%   param.numSamples   :  Kx1 or 1xK vector which gives number of samples in K subsets
%   param.w            :  weight for single shell term. 0<w<1. Default value is 0.5.
%   param.ModelFile    :  string. If set, write the model into a file 
%
%   grbParam.start     :  a warm start if provided
%   grbParam.TimeLimit :  time limit to terminate the program
%   grbParam.MIPFocus  :  MIPFocus
%
% OUTPUT
%   gradCell           :  Kx1 cell, each element is a gradient matrix.
%   grb_result         :  output from GUROBI
%
%
% Reference:
%   1. "Single- and Multiple-Shell Uniform Sampling Schemes for Diffusion MRI Using Spherical Codes", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, Peter J. Basser, IEEE Transactions on Medical Imaging, 2017.  
%   2. "Designing Single- and Multiple-Shell Sampling Schemes for Diffusion MRI Using Spherical Code", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, MICCAI 2014.  
%
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%

%%
if isfield(param,'numSamples')==0
    error('should give a number of samples!');
end
numSamples=param.numSamples;
if size(numSamples,2)==1
    numSamples=numSamples';
end
numShell = size(numSamples,2);

[tmp1, k] = size(gradAll);
if tmp1~=1
    error('wrong size of input gradients!');
end
if k~=numShell
    error('wrong size of input gradients!');
end
N = zeros(1, numShell);
for ss = 1 : numShell
    [N(ss), c] = size(gradAll{ss});
    if c~=3
        error('wrong size of input gradients!');
    end
end
if nnz(numSamples>N) > 0
    error('wrong number of samples!');
end
Nsum = sum(N);

if (isfield(param,'w')==0)
    param.w = 0.5;
end


% model
clear model 


% MIPFocus=1 gives better results than other values
if ~isfield(grbParam,'MIPFocus')
    grbParam.MIPFocus=1; 
end
if ~isfield(grbParam,'OutputFlag') || grbParam.OutputFlag
    fprintf('MIPFocus = %d\n', grbParam.MIPFocus);
end


%% MILP

% calculate the distances between points
numShellTotal = numShell+1;
minDist = zeros(1, numShellTotal);
M = zeros(1, numShellTotal);
gradAllMatrix = zeros(Nsum,3);
distAll = cell(1, numShellTotal);
tmp = cell(1, numShellTotal);
distAllM = distAll;
for ss = 1 : numShell
    distAll{ss} = real(acos(abs(gradAll{ss}*gradAll{ss}')));
    tmp{ss} = distAll{ss}(triu(distAll{ss},1)>0);
    minDist(ss) = min(tmp{ss});
    M(ss) = max(tmp{ss})-min(tmp{ss});
    
    distAllM{ss} = distAll{ss}+2*M(ss);
    
    gradAllMatrix(sum(N(1:ss-1))+1 : sum(N(1:ss)), :) = gradAll{ss};
end
distAll{numShellTotal} = real(acos(abs(gradAllMatrix*gradAllMatrix')));
tmp{numShellTotal} = distAll{numShellTotal}(triu(distAll{numShellTotal},1)>0);
tmp_aa = tmp{numShellTotal};
minDist(end) = min(tmp_aa(:));
M(end) = max(tmp_aa(:)) - min(tmp_aa(:));
distAllM{numShellTotal} = distAll{numShellTotal} + 2*M(end);

% http://mathworld.wolfram.com/SphericalCode.html
% 2*numSamples points in sphere
upperBuondEuc = sqrt(  4-csc( pi*2*numSamples./(6*(2*numSamples-2)) ).^2 );
upperBuondSpherical = acos((2.0-upperBuondEuc.*upperBuondEuc)./2.0);
upperBuondEucAll = sqrt(  4-csc( pi*2*sum(numSamples)./(6*(2*sum(numSamples)-2)) ).^2 );
upperBuondSphericalAll = acos((2.0-upperBuondEucAll.*upperBuondEucAll)./2.0);
%     upperBuondEuc = [upperBuondEuc upperBuondEucAll];
upperBuondSpherical = [upperBuondSpherical upperBuondSphericalAll];


% objective
model.obj = [param.w/numShell*ones(1,numShell), (1-param.w), zeros(1,sum(N))];
model.modelsense = 'max';

% inequality constraint
NRowEach=zeros(1, numShellTotal);
for ss = 1 : numShellTotal
    %         NRowEach(ss) = N(ss)*(N(ss)-1)/2;
    NRowEach(ss) = nnz(tmp{ss}<=upperBuondSpherical(ss));
end
%     NRowEach(end) = Nsum*(Nsum-1)/2;
rhs = zeros(sum(NRowEach),1);
A = zeros(sum(NRowEach), numShellTotal+Nsum);
A0 = zeros(sum(NRowEach), numShellTotal+Nsum);
for ss = 1 : numShellTotal
    A0(sum(NRowEach(1:ss-1))+1: sum(NRowEach(1:ss)), ss) = 1;
end


for ss = 1 : numShell
    rhsTmp = zeros(NRowEach(ss),1);
    ATmp = zeros(NRowEach(ss), N(ss));
    kk=1;
    for i = 1 : N(ss)
        for j = 1 : i-1
            if distAll{ss}(i,j) <= upperBuondSpherical(ss)
                %         fprintf('i=%d, j=%d\n', i, j);
                rhsTmp(kk) = distAllM{ss}(i,j);
                ATmp(kk,i) = M(ss);
                ATmp(kk,j) = M(ss);
                kk = kk+1;
            end
        end
    end
    rhs(sum(NRowEach(1:ss-1))+1: sum(NRowEach(1:ss))) = rhsTmp;
    A(sum(NRowEach(1:ss-1))+1: sum(NRowEach(1:ss)), (numShellTotal+1+sum(N(1:ss-1))):(numShellTotal+sum(N(1:ss))) ) = ATmp;
end
kk = sum(NRowEach(1:end-1))+1;
for i = 2 : Nsum
    for j = 1 : i-1
        if distAll{numShellTotal}(i,j) <= upperBuondSpherical(numShellTotal)
            rhs(kk) = distAllM{numShellTotal}(i,j);
            A(kk,i+numShellTotal) = M(end);
            A(kk,j+numShellTotal) = M(end);
            kk = kk+1;
        end
    end
end

A = A + A0;
% equality constraint
Aeq = zeros(numShell, numShellTotal+sum(N));
for ss = 1 : numShell
    Aeq(ss, numShellTotal+1+sum(N(1:ss-1)) : numShellTotal+sum(N(1:ss)) ) = 1;
end

model.A = [sparse(Aeq);  sparse(A)];
model.rhs = [numSamples';  rhs];
model.sense = [repmat('=',[1, numShell]), repmat('<',[1,sum(NRowEach)])];

model.vtype = [repmat('C', [1, numShellTotal]), repmat('B',[1,sum(N)])];

if (isfield(param,'lbCost')~=0)
    model.lb = [max(param.lbCost(1:numShellTotal), minDist), zeros( 1, sum(N)) ];
else
    model.lb = [minDist, zeros( 1, sum(N)) ];
end
if (isfield(param,'ubCost')~=0)
    model.ub = [min(param.ubCost, upperBuondSpherical), ones(1,sum(N)) ];
else
    model.ub = [upperBuondSpherical, ones(1,sum(N)) ];
end




%%
if isfield(param,'ModelFile') && ~strcmp(param.ModelFile,'')
    gurobi_write(model,param.ModelFile);
end
grb_result = gurobi(model,grbParam);
indexMatrix = cell(1, numShell);
for ss = 1 : numShell
    indexMatrix{ss} = grb_result.x(numShellTotal+sum(N(1:ss-1))+1 : numShellTotal+sum(N(1:ss)));
end

% grb_result.x
gradCell = cell(1, numShell);
for ss = 1 : numShell
    gradCell{ss} = gradAll{ss}(indexMatrix{ss}>1e-4,:);
end



