function performSolverPrecheck(self, checkics)  
% Perform sanity checks before solving the PDE
%

% Copyright 2015 The MathWorks, Inc.
  
% Run the checks, the PDEModel should have
%     Geometry
%     Mesh
%     BCs
%     Coefficients
%    
  if isempty(self.Geometry)
      error(message('pde:pdeModel:noGeometry'));   
  end
  
  if isempty(self.Mesh)
      error(message('pde:pdeModel:noMesh'));   
  end
  
%   if isempty(self.BoundaryConditions)
%        error(message('pde:pdeModel:noBCs'));   
%   end
  
  if isempty(self.EquationCoefficients)
       error(message('pde:pdeModel:noCoefficients'));   
  end
    
   if checkics && self.IsTimeDependent && isempty(self.InitialConditions)
        error(message('pde:pdeModel:noICs'));   
   end
  
  self.EquationCoefficients.performSolverPrechecks();
 
  if checkics && ~isempty(self.InitialConditions)
    self.InitialConditions.performSolverPrechecks();  
  end
    
end