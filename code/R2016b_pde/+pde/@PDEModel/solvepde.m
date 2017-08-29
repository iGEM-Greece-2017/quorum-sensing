function R = solvepde(self, varargin)
% solvepde - Solve the PDE
% R = solvepde(pdem) solves a stationary PDE defined by the PDEModel pdem
% and returns the result in a StationaryResults object.
%
% R = solvepde(pdem, tlist) solves a time-dependent PDE defined by the
% PDEModel pdem and returns the result in a TimeDependentResults object.
% The input argument tlist is a vector that specifies the times at which
% the solution is requested.
%
% See also pde.PDEModel, pde.PDEModel/solvepdeeig, pde.StationaryResults,
%          pde.TimeDependentResults

% Copyright 2015-2016 The MathWorks, Inc.

narginchk(1,2);
checkics = true;
performSolverPrecheck(self, checkics);
if nargin == 2 && ~self.IsTimeDependent
    warning(message('pde:pdeModel:stationaryTlistIgnored'))
elseif nargin == 1 && self.IsTimeDependent
    error(message('pde:pdeModel:timeDependentNoTlist'))
end

tlist = [];
if self.IsTimeDependent
    tlist = varargin{1};
    validateattributes(tlist,{'numeric'},{'real', 'finite', 'nonsparse', 'nonnan'});
end
u0 = [];
ut0 = [];
if self.IsTimeDependent
    if self.EquationCoefficients.mDefined()
        [u0, ut0] = self.InitialConditions.packInitialConditions();
    else
        u0 = self.InitialConditions.packInitialConditions();
    end
end
coefstruct = self.EquationCoefficients.packCoefficients();

allbcsnumeric = true;
if ~isempty(self.BoundaryConditions)
    numbcs = numel(self.BoundaryConditions.BoundaryConditionAssignments);
    for i = 1:numbcs
        thisbc = self.BoundaryConditions.BoundaryConditionAssignments(i);
        if ~thisbc.numericCoefficients()
            allbcsnumeric = false;
            break;
        end
    end
end

if self.IsTimeDependent
    [u,dudt] = self.solveTimeDependent(coefstruct, u0, ut0, tlist, ...
        self.EquationCoefficients.mDefined());
    if (self.EquationCoefficients.mDefined())
        R = pde.TimeDependentResults(self,u,tlist,dudt);
    else
        R = pde.TimeDependentResults(self,u,tlist);
    end
    
else
    if (self.EquationCoefficients.allCoefficientsNumeric() && allbcsnumeric)
        u = self.solveStationary(coefstruct);
        R = pde.StationaryResults(self,u);
    else
        [p,e,t] = self.Mesh.meshToPet();
        u0 = pdeuxpd(p,0,self.PDESystemSize);
        if ~isempty(self.InitialConditions)
            u0 = self.InitialConditions.packInitialConditions();
        end
        femodel = pde.DiscretizedPDEModel(self,p,e,t,coefstruct,u0);
        if (~femodel.IsStationaryNonlinear)
            u = self.solveStationary(coefstruct,femodel);
            R = pde.StationaryResults(self,u);
        else
            u = self.solveStationaryNonlinear(coefstruct, u0);
            R = pde.StationaryResults(self,u);
        end
    end
end

end