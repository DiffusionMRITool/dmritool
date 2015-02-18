% @MATLAB_FUNCTION_NAME@:  read gradient files into a matrix
%
% Usage:
%   gradMatrix = @MATLAB_FUNCTION_NAME@(type, number)
%   gradMatrix = @MATLAB_FUNCTION_NAME@('grad.txt')
%
% INPUT
%   type    :  gradient type or gradient file name (.txt)
%              'elec': gradients from camino based on electrostatic energy.
%              'tess': gradients based on tessellation.
%
%   number  :  if type=='elec', then number means the number of samples.
%              if type=='tess', then number means tessellation order.
%
% OUTPUT
%   grad    :  Nx3 matrix, each row is a point in sphere.
%
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%
