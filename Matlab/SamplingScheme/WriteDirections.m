function WriteDirections(grad, gradStr)
% Write gradient table to a file
%
% USAGE:
%   WriteDirections(grad, gradStr)
%
% INPUT:
%   grad     : Nx3 matrix
%   gradStr  : output file name
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%


dlmwrite(gradStr, grad, 'delimiter',' ','precision',8);

