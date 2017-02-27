function radius = CoveringRadius(grad, ind)
% Calculate the covering radius of the given samples. 
%
% USAGE:
%    radius = CoveringRadius(grad)
%    radius = CoveringRadius(grad, ind)
%
% INPUT
%    grad     :  Nx3 matrix, each row is a point in sphere.
%    ind      :  If it is given, calculate the minimal distance between 
%                the ind-th sample and the other samples in grad.
%
%
% OUTPUT
%    radius   :  covering radius, i.e., the minimal spherical distance between any two samples.
%
%
% References: 
%   1.  http://mathworld.wolfram.com/SphericalCode.html
%   2. "Designing Single- and Multiple-Shell Sampling Schemes for Diffusion MRI Using Spherical Code", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, MICCAI 2014.  
%   3. "Novel single and multiple shell uniform sampling schemes for diffusion MRI using spherical codes", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, Peter J. Basser, MICCAI 2015.  
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%

%%
if nargin == 1
    innerProductAll = abs(grad*grad');
    radius = acos(max(innerProductAll(triu(innerProductAll,1)>0)));
else
    if ind>size(grad,1) || ind<=0
        error('wrong index');
    end
    innerProduct = abs(grad(ind,:)*grad');
    radius = acos(max([innerProduct(1:ind-1), innerProduct(ind+1:end)]));
end

