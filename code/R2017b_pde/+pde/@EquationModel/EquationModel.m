classdef (Abstract) EquationModel < handle & matlab.mixin.CustomDisplay
% pde.EquationModel  An abstract base class for various PDE models.
%
%    See also CREATEPDE, pde.PDEModel, pde.ThermalModel

% Copyright 2016-2017 The MathWorks, Inc.
  
  properties
% PDESystemSize - Number of PDEs in the system
    PDESystemSize;  
  end 
 
  
  properties
% Geometry - Geometric representation of the domain
%    A scalar object that defines the geometric domain. The geometry object 
%    is created by the geometryFromEdges function for 2D or the importGeometry 
%    function for 3D. Assemblies of geometric models are currently not supported.
    Geometry;  
    
% Mesh - Discretization of the domain
%    An object that defines the meshed domain. The mesh object stores the
%    node and element data. The mesh is created using the generateMesh function.
    Mesh;              
  end
  
  properties (SetAccess='protected')
    % SolverOptions - Algorithm options for the PDE solvers 
    %
    %    See also pde.PDESolverOptions
    SolverOptions;  
  end  
      
  methods
    function obj = EquationModel(numEqns)  
      if(nargin < 1)
        numEqns = 1;
      end  
      validateattributes(numEqns,{'numeric'},...
          {'scalar','integer','positive'});
      obj.PDESystemSize = numEqns;
      obj.SolverOptions = pde.PDESolverOptions();
    end
    

    
    function set.Geometry(self, gm)
      if( ~isempty(gm) && (~isa(gm, 'pde.DiscreteGeometry') && ~isa(gm, 'pde.AnalyticGeometry')) )
        error(message('pde:pdeModel:invalidGeometry'));
      end
      if ~isscalar(gm) && ~isempty(gm)
        error(message('pde:pdeModel:noAssemSupport'));   
      end
      self.geometricConsistencyCheck(gm);

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
    
    function tf = get.IsTwoD(self)
        if isempty(self.Geometry)
            tf = [];
        else
            tf = self.Geometry.IsTwoD();
        end
    end

    % Method declaration
    gm = geometryFromEdges(self, gedges)
    gm = importGeometry(self, geofilename)
    [g,m] = geometryFromMesh(self, nodes,elems)    
    msh = generateMesh(self, varargin)  
  end
  
  properties (Dependent=true,Hidden=true)
    IsTwoD;  
  end  
  
  methods (Hidden = true)
      %Internal methods declarations
      [a,b,c] = assema(self, varargin) 
      sol = solveStationary(self, varargin)
      [sol, dsol] = solveTimeDependent(self, varargin)     
      [u, res] = solveStationaryNonlinear(self,varargin);   
      [evec, evalue] = solveEigenvalue(self,varargin);      
      femstruct = assembleFEMatricesInternal(self, varargin)
      [Q,G,H,R] = assembleBoundary(self,varargin)
      msh = generateMeshFacettedAnalytic(self,Hmax,Hmin,Hgrad,geomOrder)
      gm = analyticToDiscrete(self)
      F = facetAnalyticGeometry(self, g)
            
      function ok = geometricConsistencyCheck(~,~)
          % Derived classes can overload this method for additional checks.
          ok = true;
          return
      end
      
  end
  

  
  methods (Hidden = true, Access = protected)
     performSolverPrecheck(self,checkics)
  end
      
  methods (Hidden = true, Static = true)
     K  =  checkForSymmetry(K);
  end  
  
  methods(Hidden = true, Static  = true, Access = protected)
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
      function ok=isValidHgrad(hgrad)
        if ~isreal(hgrad) || ~isscalar(hgrad) || ischar(hgrad) || issparse(hgrad) || ~isfinite(hgrad)
            error(message('pde:pdeModel:invalidHgrad'));  
        elseif hgrad<=1 || hgrad>=2
            error(message('pde:pdeModel:invalidHgrad'));  
        end        
        ok=true;
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
      
      function ok = isValidAttributeVal(attributeVal)
          validateattributes(attributeVal,{'double','function_handle'},{})
          ok=true;
      end   
      
      function delistMesh(~,~)
%         if isvalid(self.ParentPdemodel)
%             self.ParentPdemodel.Mesh = [];
%         end
      end 
      function delistGeometry(~,~)
%         if isvalid(self.ParentPdemodel)
%             self.ParentPdemodel.Geometry = [];
%         end
      end 
  end
  
  
end

