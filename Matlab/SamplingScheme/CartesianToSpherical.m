function [theta, phi, r] = CartesianToSpherical(x,y,z)
% Convert x,y,z to theta, phi, r. 
%
% x = r.*sin(theta).*cos(phi)
% y = r.*sin(theta).*sin(phi)
% z = r.*cos(theta)
%
% USAGE:
%    [theta, phi, r] = CartesianToSpherical(x,y,z)
%
% INPUT
%    x       :  Nx1 vector with x coordinates.
%    y       :  Nx1 vector with y coordinates.
%    z       :  Nx1 vector with z coordinates.
%
% OUTPUT
%    theta, phi, r  : Nx1 vectors
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%

r = sqrt(x.^2+y.^2+z.^2);
if (r==0)
    theta = 0;
    phi = 0;
else
    theta = acos(z./r);
    phi = atan2(y,x);
end
