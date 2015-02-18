% @MATLAB_FUNCTION_NAME@:  get sh matrix based on given order and orientationMatrix
%
% Usage:
%   shMatrix = @MATLAB_FUNCTION_NAME@(order, orientationMatrix, mode)
% 
% Inputs:
%   order:               SH order
%   orientationMatrix :  orientation matrix where each row is an orientation (x,y,z)
%   mode:                'spherical' or 'cartesian' (default), for the coordinates in orientationMatrix.
%
% Outputs:
%   shMatrix:            each column is a real SH basis
%
% Copyright (c) 2013, Jian Cheng <jian.cheng.1983@gmail.com>
%


