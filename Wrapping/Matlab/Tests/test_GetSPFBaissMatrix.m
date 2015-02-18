% function test_GetSPFBaissMatrix

grad_t3 = readGrad('tess',3);
bVec = 1500*ones(size(grad_t3,1),1);

D0 = 0.7e-3;
tau = 0.02533;
param.scale = 1.0 / (8*pi^2*tau*D0);
param.original=false;
% param.original=true;

qVec = sqrt(bVec./(4*pi*pi*tau));

spfMatrix = GetSPFBasisMatrix(4,1,grad_t3,bVec,param);
