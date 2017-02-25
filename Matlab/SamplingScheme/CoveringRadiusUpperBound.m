function ub = CoveringRadiusUpperBound(N)
% calculate the upper bound () of the minimal distance for a given number of points.  
% Note: for N antipodal symmetric samples, use CoveringRadiusUpperBound(2*N).
%
% USAGE:
%    ub = CoveringRadiusUpperBound(N)
%
% INPUT
%    N  :  an integer number.
%
%
% OUTPUT
%    up :  an upper bound (spherical distance in radian)
% 
% References: 
%   1.  http://mathworld.wolfram.com/SphericalCode.html
%   2. "Designing Single- and Multiple-Shell Sampling Schemes for Diffusion MRI Using Spherical Code", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, MICCAI 2014.  
% 
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%

%%
upperBoundEuc = sqrt(  4-csc( pi*N./(6*(N-2)) ).^2 );
ub = acos((2.0-upperBoundEuc.*upperBoundEuc)./2.0);

end