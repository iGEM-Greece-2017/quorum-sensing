function performSolverPrecheck(self, checkics)
% performSolverPrecheck Perform sanity checks before solving the PDE
%

% Copyright 2015-2017 The MathWorks, Inc.

% Run the checks, the PDEModel should have
%     Geometry
%     Mesh
%     BCs
%     Coefficients
%

if isa(self,'pde.PDEModel')
    
    if isempty(self.Geometry)
        error(message('pde:pdeModel:noGeometry'));
    end
    
    if isempty(self.Mesh)
        error(message('pde:pdeModel:noMesh'));
    end
    
    
    if checkics && self.IsTimeDependent && isempty(self.InitialConditions)
        error(message('pde:pdeModel:noICs'));
    end
    
    if isempty(self.EquationCoefficients)
        error(message('pde:pdeModel:noCoefficients'));
    end
    self.EquationCoefficients.performSolverPrechecks();
    
end

if isa(self,'pde.ThermalModel')
    
    if isempty(self.Geometry)
        error(message('pde:ThermalModel:thermalModelHasNoGeom'));
    end
    
    if isempty(self.Mesh)
        error(message('pde:ThermalModel:noMesh'));
    end
    
    if isempty(self.MaterialProperties)
        error(message('pde:ThermalModel:noThermalProperties'));
    end
    self.MaterialProperties.performSolverPrechecks();
    
    if checkics && self.IsTimeDependent && isempty(self.InitialConditions)
        error(message('pde:ThermalModel:noICs'));
    end
    

end


if isa(self,'pde.StructuralModel')
    
    if isempty(self.Geometry)
        error(message('pde:StructuralModel:structuralModelHasNoGeom'));
    end
    
    if isempty(self.Mesh)
        error(message('pde:StructuralModel:noMesh'));
    end
    
    if isempty(self.MaterialProperties)
        error(message('pde:StructuralModel:noStructuralProperties'));
    end
    self.MaterialProperties.performSolverPrechecks();

    if ~isempty(self.BodyLoads)
        massEmpty = ~any([self.MaterialProperties.MaterialAssignments(:).MassDensity]);
        %gravitySpecified = any(~isempty([self.BodyLoads.BodyLoadAssignments(:).GravitationalAcceleration]));
        if massEmpty
            error(message('pde:StructuralModel:noMassDensityForGrav')); 
        end
        self.BodyLoads.performSolverPrechecks();
    end

    
end

 if checkics && ~isempty(self.InitialConditions)
      self.InitialConditions.performSolverPrechecks();
 end
    

end