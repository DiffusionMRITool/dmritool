function [grad, xopt, fopt, retcode] = OptimalSamplingSingleShellCNLO(gradInitial, param)
% update gradients from an initial gradient set, such that the updated gradients are evenly distributed. 
%
% USAGE:
%   [grad, xopt, fopt, retcode] = OptimalSamplingSingleShellCNLO(gradInitial, param)
%
% INPUT
%   gradInitial       :  N x 3 gradient matrix, where each row is a point in sphere. 
%   param.cartesian   :  true: use cartesian coordinate (with unit equality constraint);
%                        false: use spherical coordinate (no unit equality constraint)
%
% 
% OUTPUT
%   grad              :  updated N x 3 gradient matrix. 
%
%
% Reference:
%   1. "Single- and Multiple-Shell Uniform Sampling Schemes for Diffusion MRI Using Spherical Codes", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, Peter J. Basser, IEEE Transactions on Medical Imaging, 2017.  
%   2. "Novel single and multiple shell uniform sampling schemes for diffusion MRI using spherical codes", 
%       Jian Cheng, Dinggang Shen, Pew-Thian Yap, Peter J. Basser, MICCAI 2015.  
%
% Copyright (c) 2014, Jian Cheng <jian.cheng.1983@gmail.com>
%

%% sampling parameters
paramDefault.w=0.5;
paramDefault.cartesian=true;

paramDefault.maxtime=-1;
paramDefault.verbose=0;

if nargin==1
    param = paramDefault;
else
    param = CopyDefaultStructParams(paramDefault, param);
end

%%
[numSamples,c] = size(gradInitial);
if c~=3
    error('wrong size of input gradients!');
end
upperBuondEuc = sqrt(  4-csc( pi*2*numSamples./(6*(2*numSamples-2)) ).^2 );
upperBuondSpherical = acos((2.0-upperBuondEuc.*upperBuondEuc)./2.0);


%% initial
if param.cartesian
    xInit = zeros(1+3*numSamples,1);
else
    xInit = zeros(1+2*numSamples,1);
end
xInit(1) = CoveringRadius(gradInitial);
if param.cartesian
    gradAll = gradInitial';
    xInit(2:end) = gradAll(:);
else
    gradAll_sph = zeros(numSamples,2);
    [gradAll_sph(:,1), gradAll_sph(:,2)] = CartesianToSpherical(gradInitial(:,1), gradInitial(:,2), gradInitial(:,3));
    gradAll_sph = gradAll_sph';
    xInit(2:end) = gradAll_sph(:);
end


%% nlopt parameters

opt.xtol_rel = 1e-4;
opt.ftol_rel = 1e-6;

opt.maxtime = param.maxtime;

opt.algorithm = NLOPT_LD_SLSQP;


%% cost
opt.max_objective = @(x) maxmin_func(x);

if param.cartesian
    opt.h = cell(numSamples,1);
    % opt.h_tol = 1e-3*ones(numSamples,1);
    index = 1;
    for j = 1 : numSamples
        opt.h{index} = @(x) maxmin_constrains_unit(x, index);
        index = index+1;
    end
end

index = 1;
for i = 2 : numSamples
    for j = 1 : i-1
        if param.cartesian
            opt.fc{index} = @(x) maxmin_constrains_ij_cartesian(x, i,j) ;
        else
            opt.fc{index} = @(x) maxmin_constrains_ij_spherical(x, i,j) ;
        end
        index = index+1;
    end
end

% opt.fc_tol = [1e-8, 1e-8];
% opt.xtol_rel = 1e-4;

if param.cartesian
    opt.lower_bounds = [xInit(1), -1*ones(1,3*numSamples)];
    opt.upper_bounds = [upperBuondSpherical, ones(1,3*numSamples)];
else
    opt.lower_bounds = [xInit(1), repmat([0,-pi], [1,numSamples])];
    opt.upper_bounds = [upperBuondSpherical, repmat([pi/2,pi], [1,numSamples])];
end


%% solve
[xopt, fopt, retcode] = nlopt_optimize(opt, xInit);


%% output
grad = gradInitial;
if param.cartesian
    grad_i = xopt(2 : end);
    grad = reshape(grad_i, [3,numSamples])';
else
    grad_i = xopt(2 : end);
    grad_i = reshape(grad_i, [2,numSamples])';
    [grad(:,1), grad(:,2), grad(:,3)] = SphericalToCartesian(grad_i(:,1), grad_i(:,2), 1);
end
end

%% 
function [val, gradient] = maxmin_func(x)
% numTotalSamples = sum(numSamples);
val = x(1);

if (nargout > 1)
    gradient = zeros(size(x));
    gradient(1) = 1;
end
end

function [val, gradient] = maxmin_constrains_unit(x, index)
% numTotalSamples = sum(numSamples);
offset = 1;
val = x(offset+(index-1)*3+1)*x(offset+(index-1)*3+1) ...
    + x(offset+(index-1)*3+2)*x(offset+(index-1)*3+2) ...
    + x(offset+(index-1)*3+3)*x(offset+(index-1)*3+3)-1;

if (nargout > 1)
    gradient = zeros(size(x));
    gradient(offset+(index-1)*3+1) = 2*x(offset+(index-1)*3+1);
    gradient(offset+(index-1)*3+2) = 2*x(offset+(index-1)*3+2);
    gradient(offset+(index-1)*3+3) = 2*x(offset+(index-1)*3+3);
end
end

function [val, gradient] = maxmin_constrains_ij_cartesian(x, index_i, index_j)
% numTotalSamples = sum(numSamples);
offset = 1;
dotProduct = x(offset+(index_i-1)*3+1)*x(offset+(index_j-1)*3+1) ...
    + x(offset+(index_i-1)*3+2)*x(offset+(index_j-1)*3+2) ...
    + x(offset+(index_i-1)*3+3)*x(offset+(index_j-1)*3+3);
dotProduct = dotProduct/ ( norm(x(offset+(index_i-1)*3+1:offset+(index_i-1)*3+3)) * norm(x(offset+(index_j-1)*3+1:offset+(index_j-1)*3+3)));
val = x(1) - real(acos(abs(dotProduct)));

if (nargout > 1)
    gradient = zeros(size(x));
    gradient(1) = 1;
    scale = 1/sqrt(1-dotProduct*dotProduct) * sign(dotProduct);
    for ii = 1:3
        gradient(offset+(index_i-1)*3+ii) = scale * x(offset+(index_j-1)*3+ii);
        gradient(offset+(index_j-1)*3+ii) = scale * x(offset+(index_i-1)*3+ii);
    end
%     size(gradient)
end
end


function [val, gradient] = maxmin_constrains_ij_spherical(x, index_i, index_j)
% numTotalSamples = sum(numSamples);
offset = 1;

vSph_i = x(offset+(index_i-1)*2+1:offset+(index_i-1)*2+2);
vSph_j = x(offset+(index_j-1)*2+1:offset+(index_j-1)*2+2);
vCar_i = zeros(3,1); 
vCar_j = zeros(3,1);
[vCar_i(1), vCar_i(2), vCar_i(3)] = SphericalToCartesian(vSph_i(1),vSph_i(2),1);
[vCar_j(1), vCar_j(2), vCar_j(3)] = SphericalToCartesian(vSph_j(1),vSph_j(2),1);
dotProduct = vCar_i(1)*vCar_j(1) ...
    + vCar_i(2)*vCar_j(2) ...
    + vCar_i(3)*vCar_j(3);
dotProduct = dotProduct/ ( norm(vCar_i) * norm(vCar_j));
val = x(1) - real(acos(abs(dotProduct)));

if (nargout > 1)
    
    gradient = zeros(size(x));
    gradient(1) = 1;
    scale = 1/sqrt(1-dotProduct*dotProduct) * sign(dotProduct);

%     if index_i==2 && index_j==1
%         vSph_i
%         vSph_j
%         vCar_i
%         vCar_j
%         dotProduct
%         angle=acos(abs(dotProduct))*180/pi;
%         angle
%         scale
%     end

    theta_i = vSph_i(1);
    phi_i = vSph_i(2);
    theta_j = vSph_j(1);
    phi_j = vSph_j(2);
    gradient(offset+(index_i-1)*2+1) = scale * (cos(theta_i)*cos(phi_i)*vCar_j(1)+ cos(theta_i)*sin(phi_i)*vCar_j(2) - sin(theta_i)*vCar_j(3) );
    gradient(offset+(index_i-1)*2+2) = scale * (-sin(theta_i)*sin(phi_i)*vCar_j(1)+ sin(theta_i)*cos(phi_i)*vCar_j(2) );
    gradient(offset+(index_j-1)*2+1) = scale * (cos(theta_j)*cos(phi_j)*vCar_i(1)+ cos(theta_j)*sin(phi_j)*vCar_i(2) - sin(theta_j)*vCar_i(3) );
    gradient(offset+(index_j-1)*2+2) = scale * (-sin(theta_j)*sin(phi_j)*vCar_i(1)+ sin(theta_j)*cos(phi_j)*vCar_i(2) );
  
    
end
end



