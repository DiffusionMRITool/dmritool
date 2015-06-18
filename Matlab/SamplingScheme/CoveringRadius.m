function radius = CoveringRadius(grad)
% calculate the covering radius of the given samples
%
% USAGE:
%    radius = CoveringRadius(grad)
%
% INPUT
%    grad     :  Nx3 matrix, each row is a point in sphere.
%
%
% OUTPUT
%    radius   :  covering radius, i.e., the minimal spherical distance between any two samples.
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%

%%
innerProductAll = abs(grad*grad');
radius = acos(max(innerProductAll(triu(innerProductAll,1)>0)));

