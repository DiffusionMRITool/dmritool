function [gradCell, xopt, fopt, retcode] = OptimalSamplingMultiShellCNLO_singleRun(gradCellInitial, param)
% update gradients from an initial gradient set, such that the updated gradients are evenly distributed. 
%
% INPUT (required)
%   gradCellInitial      :  N_k x 3 gradient matrices, where each row is a point in sphere. k=1,2,...,K
%   param.w              :  weight for single shell term. 0<w<1. 0.5 is the default value.
% 
% INPUT (optional)
%   param.solver         :  'nlopt' (default): use sqp solver (SLSQP) in NLOPT.
%                           'matlab': use fmincon (sqp solver).
%   param.cartesian      :  true (default): use cartesian coordinate (with unit equality constraint);
%                           false: use spherical coordinate (no unit equality constraint).
% 
%   param.localCon       :  0: do not use local constraint
%                           1 (default): use local constraint for the distance between samples in results and samples in initialization to reduce the inequality constraint. 
%                                        It is faster, and with an appropriate localConAngle, it can obtain the same results as it is set as false. 
%
%   param.localConAngle  :  It is the maximal angular change between samples in results and samples in initialization. 
%                           It is used when param.localCon is true. It is automatically set based on the upper bounds, if it is not given and param.localCon=true. 
%                           If localConAngle>pi/2, then the result is the same as the results without localCon.
%
% 
% OUTPUT
%   gradCell             :  Kx1 cell, each element is a gradient matrix. 
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
paramDefault.solver='nlopt';
paramDefault.cartesian=true;
paramDefault.localCon=1;

paramDefault.maxtime=-1;
paramDefault.verbose=0;

if nargin==1
    param = paramDefault;
else
    param = CopyDefaultStructParams(paramDefault, param);
end



%%
numShells = numel(gradCellInitial);
numSamples = zeros(numShells,1);
for i = 1 : numShells
    [numSamples(i),c] = size(gradCellInitial{i});
    if c~=3
        error('wrong size of input gradients!');
    end
end

upperBoundSpherical = CoveringRadiusUpperBound(2*numSamples);
upperBoundSphericalAll = CoveringRadiusUpperBound(2*sum(numSamples));

if param.localCon
    if ~isfield(param, 'localConAngle')
        param.localConAngle = max([upperBoundSphericalAll; upperBoundSpherical]) / 1.5;
    end
end

%% initial
if param.cartesian
    xInit = zeros(numShells+1+3*sum(numSamples),1);
else
    xInit = zeros(numShells+1+2*sum(numSamples),1);
end
gradAll = zeros(sum(numSamples),3);
index=1;
for k = 1 : numShells
    xInit(k) = CoveringRadius(gradCellInitial{k});
    gradAll(index:index+numSamples(k)-1,:) = gradCellInitial{k};
    index = index+numSamples(k);
end
xInit(numShells+1) = CoveringRadius(gradAll);
if param.cartesian
    gradAll = gradAll';
    xInit(numShells+2:end) = gradAll(:);
else
    gradAll_sph = zeros(sum(numSamples),2);
    [gradAll_sph(:,1), gradAll_sph(:,2)] = CartesianToSpherical(gradAll(:,1), gradAll(:,2), gradAll(:,3));
    gradAll_sph = gradAll_sph';
    xInit(numShells+2:end) = gradAll_sph(:);
end


%% nlopt parameters

opt.xtol_rel = 1e-4;
opt.ftol_rel = 1e-6;
opt.maxtime = param.maxtime;

opt.algorithm = NLOPT_LD_SLSQP;


%% cost
opt.max_objective = @(x) maxmin_func(x, numSamples, param.w, upperBoundSpherical, upperBoundSphericalAll);

if param.cartesian
    opt.h = cell(sum(numSamples),1);
    % opt.h_tol = 1e-3*ones(sum(numSamples),1);
    index = 1;
    for k = 1 : numShells
        for j = 1 : numSamples(k)
            opt.h{index} = @(x) maxmin_constrains_unit(x, numSamples, index);
            index = index+1;
        end
    end
end

offset = numShells+1;
% i, j in the same shell (k-th shell)
index = 1;
index_shells=0;
for k = 1 : numShells
    for i = 2 : numSamples(k)
        for j = 1 : i-1
            index_i=i+index_shells;
            index_j=j+index_shells;
            % fprintf('index_i, index_j = (%d,%d)\n', index_i, index_j)
            
            if param.localCon
                dotProduct = xInit(offset+(index_i-1)*3+1)*xInit(offset+(index_j-1)*3+1) ...
                    + xInit(offset+(index_i-1)*3+2)*xInit(offset+(index_j-1)*3+2) ...
                    + xInit(offset+(index_i-1)*3+3)*xInit(offset+(index_j-1)*3+3);
                dotAbs = abs(dotProduct);
                angleTmp = upperBoundSpherical(k)+2*param.localConAngle;
                if (angleTmp>=pi/2 || dotAbs>cos(angleTmp))
                    if param.cartesian
                        opt.fc{index} = @(x) maxmin_constrains_ij_cartesian(x, numSamples, k, index_i,index_j, cos(upperBoundSpherical(k)) ) ;
                    else
                        opt.fc{index} = @(x) maxmin_constrains_ij_spherical(x, numSamples, k, index_i,index_j) ;
                    end
                    index = index+1;
                end
            else
                if param.cartesian
                    opt.fc{index} = @(x) maxmin_constrains_ij_cartesian(x, numSamples, k, index_i,index_j, cos(upperBoundSpherical(k)) ) ;
                else
                    opt.fc{index} = @(x) maxmin_constrains_ij_spherical(x, numSamples, k, index_i,index_j) ;
                end
                index = index+1;
            end
        end
    end
    index_shells = index_shells+numSamples(k);
end
% i, j in the combined shell (i, j are not in the same shell)
for k = 2 : numShells
    for k2 = 1 : k-1
        for i = 1 : numSamples(k)
            index_i = sum(numSamples(1:k-1)) + i;
            for j = 1 : numSamples(k2)
                index_j = sum(numSamples(1:k2-1)) + j;
                
                if param.localCon
                    
                    dotProduct = xInit(offset+(index_i-1)*3+1)*xInit(offset+(index_j-1)*3+1) ...
                        + xInit(offset+(index_i-1)*3+2)*xInit(offset+(index_j-1)*3+2) ...
                        + xInit(offset+(index_i-1)*3+3)*xInit(offset+(index_j-1)*3+3);
                    dotAbs = abs(dotProduct);
                    angleTmp = upperBoundSphericalAll+2*param.localConAngle;
                    if (angleTmp>=pi/2 || dotAbs>cos(angleTmp))
                        if param.cartesian
                            opt.fc{index} = @(x) maxmin_constrains_ij_cartesian(x, numSamples, numShells+1, index_i,index_j, cos(upperBoundSphericalAll) ) ;
                        else
                            opt.fc{index} = @(x) maxmin_constrains_ij_spherical(x, numSamples, numShells+1, index_i,index_j) ;
                        end
                        index = index+1;
                    end
                else
                    if param.cartesian
                        opt.fc{index} = @(x) maxmin_constrains_ij_cartesian(x, numSamples, numShells+1, index_i,index_j) ;
                    else
                        opt.fc{index} = @(x) maxmin_constrains_ij_spherical(x, numSamples, numShells+1, index_i,index_j) ;
                    end
                    index = index+1;
                    
                end
                
                
            end
        end
    end
end

% the covering radius in individual shell is no smaller than the covering radius in the combined shell
for k = 1 : numShells
    opt.fc{index} = @(x) maxmin_constrains_shell(x, k, numShells);
    index = index+1;
end

% local constraint when param.localCon is true
if param.localCon
    minDot = cos(param.localConAngle);
    for i= 1: sum(numSamples)
        if param.cartesian
            opt.fc{index} = @(x) maxmin_constrainsLocal_i_cartesian(x, xInit, numSamples, i, minDot ) ;
        else
            opt.fc{index} = @(x) maxmin_constrainsLocal_i_spherical(x, xInit, numSamples, i, minDot ) ;
        end
        index = index+1;
    end
end

% opt.fc_tol = [1e-8, 1e-8];
% opt.xtol_rel = 1e-4;

if param.cartesian
    opt.lower_bounds = [xInit(1:numShells+1)'*0.8, -1*ones(1,3*sum(numSamples))];
    opt.upper_bounds = [upperBoundSpherical', upperBoundSphericalAll,  ones(1,3*sum(numSamples))];
else
    opt.lower_bounds = [xInit(1:numShells+1)'*0.8, repmat([0,-pi], [1,sum(numSamples)])];
    opt.upper_bounds = [upperBoundSpherical', upperBoundSphericalAll, repmat([pi/2,pi], [1,sum(numSamples)])];
end

if isfield(param, 'lower_bounds')
    fprintf('lower_bounds\n');
    opt.lower_bounds(1:numShells+1)
    param.lower_bounds
    opt.lower_bounds(1:numShells+1) = param.lower_bounds;
end


%% solve using nlopt or matlab
% fprintf('start\n');
if param.verbose>0
    opt.verbose = 1;
end

if strcmp(param.solver,'nlopt')
    [xopt, fopt, retcode] = nlopt_optimize(opt, xInit);
elseif strcmp(param.solver,'matlab')
    % solve matlab
    nonlcon = @(x) nonlcon_matlab(x, opt );
    fun = @(x)(-opt.max_objective(x) );
    if param.verbose>0
        opt_fmincon = optimoptions('fmincon','Display','iter', 'Algorithm','sqp');
    else
        opt_fmincon = optimoptions('fmincon','Display','final','Algorithm','sqp');
    end
    [xopt, fopt, retcode, outinfo] = fmincon(fun,xInit,[],[],[],[],opt.lower_bounds,opt.upper_bounds,nonlcon,opt_fmincon);
end
% fprintf('end\n');

%% output
gradCell = gradCellInitial;
index_shells=0;
for k = 1 : numShells
    if param.cartesian
        grad_i = xopt(numShells+1+index_shells+1 : numShells+1+index_shells + 3*numSamples(k));
        gradCell{k} = reshape(grad_i, [3,numSamples(k)])';
        index_shells = index_shells + 3*numSamples(k);
    else
        grad_i = xopt(numShells+1+index_shells+1 : numShells+1+index_shells + 2*numSamples(k));
        grad_i = reshape(grad_i, [2,numSamples(k)])';
        [gradCell{k}(:,1), gradCell{k}(:,2), gradCell{k}(:,3)] = SphericalToCartesian(grad_i(:,1), grad_i(:,2), 1);
        index_shells = index_shells + 2*numSamples(k);
    end
end
end

%% 
function [val, gradient] = maxmin_func(x, numSamples, w, upperBoundSpherical, upperBoundSphericalAll)
% numTotalSamples = sum(numSamples);
numShells = numel(numSamples);
val = 0;
for i = 1:numShells
%     val = val + w/(numShells*upperBoundSpherical(i))*x(i);
    val = val + w/(numShells)*x(i);
end
% val = val + (1-w)/upperBoundSphericalAll*x(numShells+1);
val = val + (1-w)*x(numShells+1);

if (nargout > 1)
    gradient = zeros(size(x));
    for i = 1:numShells
%         gradient(i) = w/(numShells*upperBoundSpherical(i));
        gradient(i) = w/numShells;
    end
%     gradient(numShells+1)=(1-w)/upperBoundSphericalAll;
    gradient(numShells+1)=(1-w);
end
end

function [val, gradient] = maxmin_constrains_unit(x, numSamples, index)
% numTotalSamples = sum(numSamples);
numShells = numel(numSamples);
offset = numShells+1;
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

function [val, gradient] = maxmin_constrains_ij_cartesian(x, numSamples, index_shell, index_i, index_j, dotMinBound)
% numTotalSamples = sum(numSamples);
numShells = numel(numSamples);
offset = numShells+1;
dotProduct = x(offset+(index_i-1)*3+1)*x(offset+(index_j-1)*3+1) ...
    + x(offset+(index_i-1)*3+2)*x(offset+(index_j-1)*3+2) ...
    + x(offset+(index_i-1)*3+3)*x(offset+(index_j-1)*3+3);
dotProduct = dotProduct/ ( norm(x(offset+(index_i-1)*3+1:offset+(index_i-1)*3+3)) * norm(x(offset+(index_j-1)*3+1:offset+(index_j-1)*3+3)));
dotAbs = abs(dotProduct);

if nargin<6
    dotMinBound = 0;
end

if dotAbs > dotMinBound
    val = x(index_shell) - real(acos(dotAbs));
    % val = dotProduct*dotProduct - cos(x(index_shell))*cos(x(index_shell));
else
    val = -1;
end


if nargout > 1
    gradient = zeros(size(x));
    
    if dotAbs > dotMinBound
        
        gradient(index_shell) = 1;
        scale = 1/sqrt(1-dotProduct*dotProduct) * sign(dotProduct);
        for ii = 1:3
            gradient(offset+(index_i-1)*3+ii) = scale * x(offset+(index_j-1)*3+ii);
            gradient(offset+(index_j-1)*3+ii) = scale * x(offset+(index_i-1)*3+ii);
        end
        
        %     gradient(index_shell) = 2*cos(x(index_shell))*sin(x(index_shell));
        %     scale = 2*dotProduct;
        %     for ii = 1:3
        %         gradient(offset+(index_i-1)*3+ii) = scale * x(offset+(index_j-1)*3+ii);
        %         gradient(offset+(index_j-1)*3+ii) = scale * x(offset+(index_i-1)*3+ii);
        %     end
    end

end
end


function [val, gradient] = maxmin_constrains_ij_spherical(x, numSamples, index_shell, index_i, index_j)
% numTotalSamples = sum(numSamples);
numShells = numel(numSamples);
offset = numShells+1;

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
val = x(index_shell) - real(acos(abs(dotProduct)));

if (nargout > 1)
    
    gradient = zeros(size(x));
    gradient(index_shell) = 1;
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

function [val, gradient] = maxmin_constrains_shell(x, index_shell, numShells)
val = x(numShells+1) - x(index_shell);
if (nargout > 1)
    gradient = zeros(size(x));
    gradient(numShells+1) = 1;
    gradient(index_shell) = -1;
end
end

function [val, gradient] = maxmin_constrainsLocal_i_cartesian(x, xInit, numSamples, index_i, dotMinBound)
% numTotalSamples = sum(numSamples);
numShells = numel(numSamples);
offset = numShells+1;
dotProduct = x(offset+(index_i-1)*3+1)*xInit(offset+(index_i-1)*3+1) ...
    + x(offset+(index_i-1)*3+2)*xInit(offset+(index_i-1)*3+2) ...
    + x(offset+(index_i-1)*3+3)*xInit(offset+(index_i-1)*3+3);
dotProduct = dotProduct/ ( norm(x(offset+(index_i-1)*3+1:offset+(index_i-1)*3+3)) * norm(xInit(offset+(index_i-1)*3+1:offset+(index_i-1)*3+3)));
% dotAbs = abs(dotProduct);

dotMin2 = dotMinBound*dotMinBound;

val = dotMin2-dotProduct*dotProduct;

if nargout > 1
    gradient = zeros(size(x));
    
    scale = -2*dotProduct;
    for ii = 1:3
        gradient(offset+(index_i-1)*3+ii) = scale * xInit(offset+(index_i-1)*3+ii);
    end
    
end
end

function [c, ceq] = nonlcon_matlab(x, opt)
    c = zeros(length(opt.fc),1);
    ceq = zeros(length(opt.h),1);
    
    for i = 1: length(opt.fc)
        c(i) = opt.fc{i}(x);
    end
    for i = 1: length(opt.h)
        ceq(i) = opt.h{i}(x);
    end
end
