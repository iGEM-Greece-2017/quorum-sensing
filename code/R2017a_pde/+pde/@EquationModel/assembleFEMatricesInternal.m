function FEMatricesOut = assembleFEMatricesInternal(self, varargin)
% assembleFEMatricesInternal - Internal function to assemble the intermediate FE matrices

% Copyright 2015-2016 The MathWorks, Inc.

narginchk(1, 2)

BCEnforcementOption = 'none';
if nargin == 2
    BCEnforcementOption = varargin{1};
end

performSolverPrecheck(self, false);
[p,e,t] = self.Mesh.meshToPet();

coefstruct = self.EquationCoefficients.packCoefficients();
u0 = pdeuxpd(p,0,self.PDESystemSize);
tdummy = 0:0.1:1;
femodel = pde.DynamicDiscretizedPDEModel(self,p,e,t,coefstruct,u0,tdummy,self.EquationCoefficients.mDefined());

if (femodel.IsSpatialCoefficientsNonlinear || femodel.vd)
    error(message('pde:pdeModel:coefsDependsOnUorT'))
end

FEMatrices.K = femodel.K;
FEMatrices.A = femodel.A;
FEMatrices.F = femodel.F;
FEMatrices.Q = femodel.Q;
FEMatrices.G = femodel.G;
FEMatrices.H = femodel.H;
FEMatrices.R = femodel.R;
FEMatrices.M = femodel.Mass;


switch BCEnforcementOption
    case 'nullspace'
        [FEMatricesOut.Kc,FEMatricesOut.Fc,FEMatricesOut.B,FEMatricesOut.ud, FEMatricesOut.M] ...
            = applyBCNullSpace(FEMatrices);
    case 'stiff-spring'
        [FEMatricesOut.Ks,FEMatricesOut.Fs,FEMatricesOut.M] = applyBCStiffSpring(FEMatrices,p);
    otherwise
        FEMatricesOut = FEMatrices;
end

end

function [Kc, Fc, null, ud, M] = applyBCNullSpace(FEMatrices)
K = FEMatrices.K;
A = FEMatrices.A;
F = FEMatrices.F;
Q = FEMatrices.Q;
G = FEMatrices.G;
H = FEMatrices.H;
R = FEMatrices.R;
M = FEMatrices.M;

[null,orth]=pdenullorth(H);
if size(orth,2)==0
    ud=zeros(size(K,2),1);
else
    ud=full(orth*((H*orth)\R));
end
Kc=K+A+Q;
Fc=null'*((F+G)-Kc*ud);
Kc=null'*Kc*null;
Kc = pde.PDEModel.checkForSymmetry(Kc);
M = null'*M*null;

end

function [Ks, Fs, M] = applyBCStiffSpring(FEMatrices,p)

K = FEMatrices.K;
A = FEMatrices.A;
F = FEMatrices.F;
Q = FEMatrices.Q;
G = FEMatrices.G;
H = FEMatrices.H;
R = FEMatrices.R;
M = FEMatrices.M;

if nnz(H)==0 % No constraints
    K=K+A+Q;
    Ks = pde.PDEModel.checkForSymmetry(K);
    Fs=F+G;
    return
end

% Scaling
d1=full(diag(H*H'));
l=find(d1==0);
d1(l)=Inf*ones(size(l));
id1=diag(sparse(1./d1));
% Form normal equations
R=H'*id1*R;
H=H'*id1*H;

% Incorporate constraints (Dirichlet conditions)
% Find a positive spring constant of suitable size.

np=size(p,2);
tscale=1/sqrt(np);
kscale=max(max(abs(K)));
mscale=max(max(abs(A)));
qscale=max(max(abs(Q)));
hscale=1;

% Spring constant
spc=max(max(kscale*tscale^2,mscale),qscale*tscale)/(hscale*tscale^3);

K=K+A+Q+spc*H;
Ks =pde.PDEModel.checkForSymmetry(K);
Fs= F+G+spc*R;
M = M +spc*H;
M =pde.PDEModel.checkForSymmetry(M);

end


