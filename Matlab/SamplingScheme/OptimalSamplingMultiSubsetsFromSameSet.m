function [gradCell,grb_result,indexMatrix] = OptimalSamplingMultiSubsetsFromSameSet(gradAll, param, grbParam)
% get K subsets of gradients from one given gradient set, such that
% gradients in each subset are evenly distributed.
%
% USAGE:
%   [gradCell,grb_result,indexMatrix] = OptimalSamplingMultiSubsetsFromSameSet(gradAll, param, grbParam)
%
% INPUT
%   gradAll           :  Nx3 gradient matrix, where each row is a point in sphere. 
%   param.numSamples  :  Kx1 or 1xK vector which gives number of samples in K subsets
%   param.w           :  weight for single shell term. 0<w<1
%   param.sos         :  1 (default), use sos constraint
%                        0, do not use sos
%
%   grbParam.start     :  a warm start if provided
%   grbParam.TimeLimit :  time limit to terminate the program
%   grbParam.MIPFocus  :  MIPFocus
% 
% OUTPUT
%   gradCell          :  Kx1 cell, each element is a gradient matrix. 
%   grb_result        :  output from GUROBI
%
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%

%%
[N, c] = size(gradAll);
if c~=3
    error('wrong size of input gradients!');
end

if isfield(param,'numSamples')==0
    error('should give a number of samples!');
end
numSamples=param.numSamples;
if size(numSamples,2)==1
    numSamples=numSamples';
end
if sum(numSamples)>N
    error('wrong number of samples!');
end
numShell = size(numSamples,2);

if (isfield(param,'w')==0)
    param.w = 0.5;
end

if (isfield(param,'sos')==0)
    param.sos = 0;
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
distAll = real(acos(abs(gradAll*gradAll')));

tmp = distAll(triu(distAll,1)>0);
M = max(tmp)-min(tmp);

distAllM = distAll+2*M;

% http://mathworld.wolfram.com/SphericalCode.html
% 2*numSamples points in sphere
upperBuondEuc = sqrt(  4-csc( pi*2*numSamples./(6*(2*numSamples-2)) ).^2 );
upperBuondSpherical = acos((2.0-upperBuondEuc.*upperBuondEuc)./2.0);
upperBuondEucAll = sqrt(  4-csc( pi*2*sum(numSamples)./(6*(2*sum(numSamples)-2)) ).^2 );
upperBuondSphericalAll = acos((2.0-upperBuondEucAll.*upperBuondEucAll)./2.0);
upperBuondSpherical = [upperBuondSpherical upperBuondSphericalAll];

% objective
model.obj = [param.w/numShell*ones(1,numShell), (1-param.w), zeros(1,numShell*N)];
model.modelsense = 'max';

% inequality constraint
numShellTotal = numShell+1;
NRowEach=zeros(1, numShellTotal);
for ss = 1 : numShellTotal
    NRowEach(ss) = nnz(tmp<=upperBuondSpherical(ss));
end
A = zeros(sum(NRowEach), numShellTotal+numShell*N);
A0 = zeros(sum(NRowEach), numShellTotal+numShell*N);
rhs = zeros(sum(NRowEach),1);
for ss = 1 : numShellTotal
    A0(sum(NRowEach(1:ss-1))+1: sum(NRowEach(1:ss)), ss) = 1;
end

for ss = 1 : numShell
    rhsTmp = zeros(NRowEach(ss),1);
    ATmp = zeros(NRowEach(ss), N);
    kk=1;
    for i = 2 : N
        for j = 1 : i-1
            if distAll(i,j) <= upperBuondSpherical(ss)
                %         fprintf('i=%d, j=%d\n', i, j);
                rhsTmp(kk) = distAllM(i,j);
                ATmp(kk,i) = M;
                ATmp(kk,j) = M;
                kk = kk+1;
            end
        end
    end
    rhs(sum(NRowEach(1:ss-1))+1: sum(NRowEach(1:ss))) = rhsTmp;
    A(sum(NRowEach(1:ss-1))+1: sum(NRowEach(1:ss)), (numShellTotal+1+(ss-1)*N):(numShellTotal+ss*N) ) = ATmp;
end
kk = sum(NRowEach(1:numShell))+1;
for i = 2 : N
    for j = 1 : i-1
        if distAll(i,j) <= upperBuondSpherical(numShell+1)
            rhs(kk) = distAllM(i,j);
            % because of sos constraint, we can set all i and j across shells
            A(kk,(i+numShellTotal):N:end) = M;
            A(kk,(j+numShellTotal):N:end) = M;
            kk = kk+1;
        end
    end
end

A = A + A0;
% equality constraint
Aeq = zeros(numShell, numShellTotal+numShell*N);
for ss = 1 : numShell
    Aeq(ss, numShellTotal+1+(ss-1)*N : numShellTotal+ss*N ) = 1;
end

if param.sos==0
    % consistency constraint
    Acons = [zeros(N, numShellTotal), repmat(eye(N), [1,numShell]) ];
    
    %     % relation between minimal angles in shells
    %     Aminangle = [-1*eye(numShell), ones(numShell,1), zeros(numShell, numShell*N)];
    %
    %     model.A = [sparse(Aeq); sparse(Acons); Aminangle; A];
    %     model.rhs = [numSamples'; ones(N,1); zeros(numShell,1); rhs];
    %     model.sense = [repmat('=',[1, numShell]), repmat('<',[1,N]), repmat('<',[1,numShell]), repmat('<',[1,numShellTotal*NRowEach])];
    
    model.A = [sparse(Aeq); sparse(Acons); sparse(A)];
    model.rhs = [numSamples'; ones(N,1); rhs];
    model.sense = [repmat('=',[1, numShell]), repmat('<',[1,N]), repmat('<',[1,sum(NRowEach)])];
else
    model.A = [sparse(Aeq); sparse(A)];
    model.rhs = [numSamples'; rhs];
    model.sense = [repmat('=',[1, numShell]), repmat('<',[1,sum(NRowEach)])];
    for nn = 1 : N
        model.sos(nn).type = 1;
        model.sos(nn).index = numShellTotal+nn : N : numShellTotal+N*numShell;
    end
end

model.vtype = [repmat('C', [1, numShellTotal]), repmat('B',[1,numShell*N])];

if (isfield(param,'lbCost')~=0)
    model.lb = [max(param.lbCost, min(tmp)), zeros( 1, numShell*N) ];
else
    model.lb = [repmat(min(tmp), [1, numShellTotal]), zeros( 1, numShell*N) ];
end
if (isfield(param,'ubCost')~=0)
    model.ub = [min(param.ubCost, upperBuondSpherical), ones(1,numShell*N) ];
else
    model.ub = [upperBuondSpherical, ones(1,numShell*N) ];
end




%%
grb_result = gurobi(model, grbParam);
indexMatrix = reshape(grb_result.x(numShellTotal+1:end), [N, numShell] );

% grb_result.x
gradCell = cell(numShell);
for ss = 1 : numShell
    gradCell{ss} = gradAll(indexMatrix(:,ss)>1e-4,:);
end



