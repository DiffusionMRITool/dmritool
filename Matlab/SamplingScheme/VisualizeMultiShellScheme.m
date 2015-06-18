function VisualizeMultiShellScheme(varargin)
% visualize single or multiple shell schemes
%
% USAGE:
%   VisualizeMultiShellScheme(grad1)
%   VisualizeMultiShellScheme(grad1, grad2, grad3)
%
% INPUT
%   varargin :  One or more matrix for schemes in different shells. 
%               Each matrix has size N_ix3, N_i is the number of samples in the i-th shell.
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%


colors{1} = 'ro';
colors{2} = 'go';
colors{3} = 'bo';
colors{4} = 'co';
colors{5} = 'bo';
ss=100;

figure; 
axis equal; 
allDirection = [];
for i = 1 : nargin
    direction = varargin{i};
    scatter3([direction(:,1);-direction(:,1)], [direction(:,2); -direction(:,2)], [direction(:,3); -direction(:,3)],ss, colors{i}, 'filled'); hold on;
    allDirection=[allDirection; direction; -direction; ];
end

DT=DelaunayTri(allDirection);
tetramesh(DT, 'FaceColor', [0.8 0.8 0.8]); axis equal; 
rotate3d

