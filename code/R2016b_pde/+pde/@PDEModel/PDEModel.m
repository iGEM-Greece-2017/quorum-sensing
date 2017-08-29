classdef (Sealed) PDEModel < handle
% pde.PDEMODEL  A container that represents a PDE model
%    A PDEModel is a container that stores the geometry, mesh, and boundary
%    conditions that are required to define a PDE. A PDEModel is generally 
%    not created directly, it is created using the CREATEPDE function. 
%    The PDEModel provides functions that allow you to: define 2D geometry 
%    or import 3D geometry, generate and store a mesh for these geometries, 
%    and apply and store boundary conditions that define the behavior at 
%    the boundary of the geometry. The properties of the PDEModel allow you 
%    to access this information easily.
%
%    PDEM = pde.PDEMODEL(PDESystemSize) creates a PDE model for a system of
%    PDESystemSize partial differential equations. If PDESystemSize is
%    omitted, it defaults to one (scalar system).
%
% PDEModel methods:
%    geometryFromEdges       - Create 2D geometry from DECSG or PDEGEOM
%    importGeometry          - Create 3D geometry from file (STL only)
%    geometryFromMesh        - Create 3D geometry from a triangulated mesh
%    specifyCoefficients     - Specify all PDE coefficients over a domain or subdomain
%    applyBoundaryCondition  - Apply a boundary condition to the geometry
%    setInitialConditions    - Set the initial conditions for the PDE
%    generateMesh            - Generate a mesh from the geometry
%    solvepde                - Solve the PDE
%    solvepdeeig             - Solve the PDE eigenvalue problem
%
% PDEModel properties:
%    PDESystemSize        - Number of PDEs in the system
%    IsTimeDependent      - True if the PDE is time-dependent
%    Geometry             - Geometric representation of the domain
%    EquationCoefficients - Equation coefficients within the domain
%    BoundaryConditions   - Boundary conditions applied to the geometry
%    InitialConditions    - Initial conditions/guess
%    Mesh                 - Discretization of the domain
%    SolverOptions        - Algorithm options for the PDE solvers
%
%
%    See also CREATEPDE, pde.GeometricModel, pde.BoundaryCondition, pde.FEMesh

% Copyright 2014-2016 The MathWorks, Inc.
  
  properties
      
% PDESystemSize - Number of PDEs in the system
    PDESystemSize;  
        
  end 
 
  properties (Dependent = true, SetAccess='private')
% IsTimeDependent - True if the PDE is time-dependent
    IsTimeDependent;  
  end
  
  properties
% Geometry - Geometric representation of the domain
%    A scalar object that defines the geometric domain. The geometry object 
%    is created by the geometryFromEdges function for 2D or the importGeometry 
%    function for 3D. Assemblies of geometric models are currently not supported.
    Geometry;  
    
% EquationCoefficients - Equation coefficients within the domain
%    An object that contains the equation coefficient assignments within the 
%    geometric domain. Equation coefficients are specified and assigned to 
%    the geometry using the specifyCoefficients function. If multiple assignments
%    are made, they are all recorded within this object.
%
%    See also pde.CoefficientAssignmentRecords/findCoefficients,
%             pde.CoefficientAssignmentRecords/globalCoefficients
    EquationCoefficients;
    
% BoundaryConditions - Boundary conditions applied to the geometry
%    An object that contains the boundary condition assignments on the
%    boundary of the geometric. Boundary conditions are specified and assigned to 
%    the geometry using the applyBoundaryCondition function. If multiple assignments
%    are made, they are all recorded within this object.
%
%    See also pde.BoundaryConditionRecords/findBoundaryConditions
    BoundaryConditions;

% InitialConditions - Initial conditions/guess
%    An object that contains the initial condition assignments within the 
%    geometric domain. Initial conditions are specified and assigned to 
%    the geometry using the setInitialConditions function. If multiple assignments
%    are made, they are all recorded within this object.
%
%    See also pde.InitialConditionsRecords/findInitialConditions,
%             pde.GeometricInitialConditions, pde.NodalInitialConditions
    InitialConditions;     
    
% Mesh - Discretization of the domain
%    An object that defines the meshed domain. The mesh object stores the
%    node and element data. The mesh is created using the generateMesh function.
    Mesh;              
  end
  
  properties (SetAccess='private')
    % SolverOptions - Algorithm options for the PDE solvers 
    %
    %    See also pde.PDESolverOptions
    SolverOptions;  
  end  
      
  methods
    function obj = PDEModel(numEqns)  
      if(nargin < 1)
        numEqns = 1;
      end  
      validateattributes(numEqns,{'numeric'},...
          {'scalar','integer','positive'});
      obj.PDESystemSize = numEqns;
      obj.IsTwoD = true;
      obj.SolverOptions = pde.PDESolverOptions();
    end
    
    function set.BoundaryConditions(self, bc)
      if(~isempty(bc) && ~isa(bc, 'pde.BoundaryConditionRecords'))
        error(message('pde:pdeModel:invalidBoundaryConditions'));
      end
      if ~isscalar(bc) && ~isempty(bc)
          bcTemp = bc(1);
          for i = 2:length(bc)
              bcTemp.BoundaryConditionAssignments(end+1) =  bc(i).BoundaryConditionAssignments;
          end
          bc = bcTemp;
      end
      if isa(bc,'pde.BoundaryConditionRecords') && isempty(bc.ParentPdemodel.Geometry)
        bc.ParentPdemodel = self;
      end
      self.BoundaryConditions = bc;
    end
    
    function set.Geometry(self, gm)
      if( ~isempty(gm) && (~isa(gm, 'pde.DiscreteGeometry') && ~isa(gm, 'pde.AnalyticGeometry')) )
        error(message('pde:pdeModel:invalidGeometry'));
      end
      if ~isscalar(gm) && ~isempty(gm)
        error(message('pde:pdeModel:noAssemSupport'));   
      end
      if isa(gm, 'pde.DiscreteGeometry')
          self.IsTwoD = false; %#ok
      else
          self.IsTwoD = true; %#ok
      end          
      self.Geometry = gm;
    end
    
    function set.Mesh(self, msh)
      if(~isempty(msh) && ~isa(msh, 'pde.FEMesh'))
        error(message('pde:pdeModel:invalidMesh'));
      end
      if ~isscalar(msh) && ~isempty(msh)
        error(message('pde:pdeModel:noAssemSupport'));   
      end
      self.Mesh = msh;
    end        
    
    function set.PDESystemSize(self, sz)
      validateattributes(sz,{'numeric'},{'scalar','integer','positive'});
      self.PDESystemSize = sz;
    end
  
    function set.EquationCoefficients(self, eqc)
      if(~isempty(eqc) && ~isa(eqc, 'pde.CoefficientAssignmentRecords'))
        error(message('pde:pdeModel:invalidEquationCoefficients'));
      end
      if ~isscalar(eqc) && ~isempty(eqc)
        error(message('pde:pdeModel:nonScalarEquationCoefficients'));   
      end
      self.EquationCoefficients = eqc;
    end       
    
    
    function set.InitialConditions(self, icr)
      if(~isempty(icr) && ~isa(icr, 'pde.InitialConditionsRecords'))
        error(message('pde:pdeModel:invalidInitialConditions'));
      end
      if ~isscalar(icr) && ~isempty(icr)
        error(message('pde:pdeModel:nonScalarInitialConditions'));   
      end
      self.InitialConditions = icr;
    end       
    
     function tf = get.IsTimeDependent(self)
        tf = false;
        if ~isempty(self.EquationCoefficients)  
            tf = self.EquationCoefficients.timeDependent();
        end            
     end
    % Method declaration
         
    gm = geometryFromEdges(self, gedges)
    gm = importGeometry(self, geofilename)
    [g,m] = geometryFromMesh(self, nodes,elems)    
    coef = specifyCoefficients(self, varargin)
    bc = applyBoundaryCondition(self,varargin)     
    ic = setInitialConditions(self, varargin)
    msh = generateMesh(self, varargin)  
    sol = solvepde(self, varargin)
    [eval, evec] = solvepdeeig(self, varargin)             
  end
  
  properties (Hidden = true, SetAccess='private')
    IsTwoD;  
  end  
  
  methods (Hidden = true)
      [a,b,c] = assema(self, varargin) 
      sol = solveStationary(self, varargin)
      [sol, dsol] = solveTimeDependent(self, varargin)     
      [u, res] = solveStationaryNonlinear(self,varargin);   
      [evec, evalue] = solveEigenvalue(self,varargin);      
      femstruct = assembleFEMatricesInternal(self, varargin)
      [Q,G,H,R] = assembleBoundary(self,varargin)
  end  
  
  methods (Hidden = true, Access = private)
     performSolverPrecheck(self,checkics)
  end
      
  methods (Hidden = true, Static = true)
     K  =  checkForSymmetry(K);
  end  
  
  methods(Static, Access = private)
      function ok = isValidHmax(hval)
        if ~isreal(hval) || ~isscalar(hval) || ischar(hval) || hval < 0 || issparse(hval) || ~isfinite(hval)
            error(message('pde:pdeModel:invalidHmax'));   
        end     
        ok = true;
      end
      function ok = isValidHmin(hval)
        if ~isreal(hval) || ~isscalar(hval) || ischar(hval) || hval < 0 || issparse(hval) || ~isfinite(hval)
            error(message('pde:pdeModel:invalidHmin'));   
        end     
        ok = true;
      end
      function ok = isValidGeomOrder(go)
        if(~ischar(go))
          error(message('pde:pdeModel:invalidGeomOrder'))
        end
        nc = numel(go);        
        if ~(strncmpi(go,'linear',nc) || strncmpi(go,'quadratic',nc))
          error(message('pde:pdeModel:invalidGeomOrder'))
        end
        ok=true;
      end
  end
  
 
  methods(Static, Hidden = true)
      function ok = isValidEntityID(entityids) 
          validateattributes(entityids,{'numeric'},{'integer','positive', 'nonzero','real', 'nonsparse'});             
          ok=true;
      end   
      
      function delistMesh(self,~)
%         if isvalid(self.ParentPdemodel)
%             self.ParentPdemodel.Mesh = [];
%         end
      end 
      function delistGeometry(self,~)
%         if isvalid(self.ParentPdemodel)
%             self.ParentPdemodel.Geometry = [];
%         end
      end 
  end
  methods(Hidden = true, Access = {?pde.CoefficientAssignmentRecords})
    function delistCoefficientAssignments(self)
        if isvalid(self)
            self.EquationCoefficients = [];        
        end
    end          
  end
  
  methods(Hidden = true, Access = {?pde.InitialConditionsRecords})
    function delistInitialConditions(self)
         if isvalid(self)
            self.InitialConditions = [];  
         end
    end         
  end
  
  methods(Hidden = true, Access = {?pde.BoundaryConditionRecords})
    function delistBoundaryConditions(self)
         if isvalid(self)
            self.BoundaryConditions = [];  
         end
    end         
  end
  
  
end

