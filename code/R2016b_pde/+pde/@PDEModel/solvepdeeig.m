function res = solvepdeeig(self, evalinterval)
% solvepdeeig - Solve the PDE eigenvalue problem
% R = solvepdeeig(pdem, EVI) solves a PDE eigenvalue problem defined by the 
% PDEModel pdem and returns the result in an EigenResults object. The input
% EVI is a two-element vector that defines the range of eigenvalues of
% interest, the lower bound of this range may be set to -Inf. solvepde 
% returns all eigenvalues within the specified range.
%
% See also pde.PDEModel, pde.PDEModel/solvepde

% Copyright 2015 The MathWorks, Inc.

validateattributes(evalinterval,{'numeric'},{'real', 'nonsparse',...
                             'vector','numel',2,'nondecreasing'},2);

performSolverPrecheck(self, false);

if (self.EquationCoefficients.mDefined() && self.EquationCoefficients.dDefined())
    error(message('pde:pdeModel:NoEigMandD'))
end

if ~(self.EquationCoefficients.mDefined() || self.EquationCoefficients.dDefined())
    error(message('pde:pdeModel:NeedMorDForEig'))
end

coefstruct = self.EquationCoefficients.packCoefficients();

% Ignore forcing when solving eigenvalue problem.
% Assign appropriate sized zeros to f-coefficient so that global F is not
% computed in calls to built-in functions.
for sd = 1:length(coefstruct.Coefficients.c) %Assign zeros for all subdomains
    coefstruct.Coefficients.f{sd} = zeros(self.PDESystemSize,1);
end

% Compute BC matrices.
[p,e,t] = self.Mesh.meshToPet();

u0 = pdeuxpd(p,0,self.PDESystemSize);
tdummy = 0:0.1:1;
femodel = pde.DynamicDiscretizedPDEModel(self,p,e,t,coefstruct,u0,tdummy,self.EquationCoefficients.mDefined());

if (femodel.IsSpatialCoefficientsNonlinear || femodel.vd)
    error(message('pde:pdeModel:coefsNonlinear'))
end

K = femodel.K;
A = femodel.A;
Q = femodel.Q;
H = femodel.H;
M = femodel.Mass;

% Impose BC.
[null,~]=pdenullorth(H);
KK=K+Q+A;
K=null'*KK*null;
B = null;
M = B'*M*B;

if (issymmetric(K)&&issymmetric(M)&&evalinterval(1)~=-Inf)
    % Positive definite is checked inside sptarn
    spd=1;
else
    spd=0;
end

[v,l,ires]=sptarn(K,M,evalinterval(1),evalinterval(2),spd, 'Stats', 'on');

if ires<0,
    error(message('pde:pdeModel:MoreEigenvaluesMayExist'));
end

if ~isempty(v)
    v=B*v;
end
res = pde.EigenResults(self,v,l);

end
