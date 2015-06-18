function grad = ReadDirections(filename)
% read direction file into a matrix
%
% USAGE:
%   grad = ReadDirections(filename) 
%
% INPUT
%   filename :   direction file name (.txt)
%
%
% OUTPUT
%   grad     :  Nx3 matrix, each row is a point in sphere.
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%


%%
fid = fopen(filename);
gradC = textscan(fid, '%f %f %f');
fclose(fid);

% for some directions which the first line shows the number of directions
grad = [gradC{1} gradC{2} gradC{3}];
if int16(grad(1,1)) == grad(1,1) && grad(1,2)~=grad(1,2)
    grad = grad(2:end, :);
end

normFactor = sqrt(sum(grad.^2,2));
for i = 1 : size(grad,1)
    if (normFactor(i)>0)
        grad(i,:) = grad(i,:) / normFactor(i);
    end
end

