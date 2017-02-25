function gradAll = GetGradAllShells(gradCell)
% Get Nx3 gradient table from gradients in all shells. 
%
% USAGE:
%    gradAll = GetGradAllShells(gradCell)
%
% INPUT
%    gradCell  :  a cell array with M gradient matrix from M shells.
%
%
% OUTPUT
%    gradAll   :  Nx3 gradient matrix from all M shells.
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%

numShells = numel(gradCell);
numSamples = zeros(numShells,1);
for i = 1 : numShells
    [numSamples(i),~] = size(gradCell{i});
end

gradAll = zeros(sum(numSamples),3);
index=1;
for k = 1 : numShells
    gradAll(index:index+numSamples(k)-1,:) = gradCell{k};
    index = index+numSamples(k);
end

end
