function ic = setInitialConditions(self, varargin) 
%setInitialConditions - Set initial conditions or initial guess
%   Set the initial conditions (ICs) for time dependent problems or the
%   initial guess for nonlinear stationary problems. The PDE toolbox can 
%   solve equations of the form:
%    
%              m*d^2u/dt^2 + d*du/dt - div(c*grad(u)) + a*u = f
%  
%   The equation to solve is defined in terms of the coefficients m, d, c, 
%   a, f. If the PDE to solve has an 'm' or 'd' coefficient defined, then
%   the PDE is time-dependent and Initial Conditions must be specified. 
%   If both the 'm' and 'd' coefficients are numeric, scalar, and zero and 
%   the 'c', 'a', or 'f' coefficients are a function of u or grad(u), then 
%   the PDE is a nonlinear stationary problem. An optional initial guess can 
%   be specified, otherwise a default initial guess of zero is assumed.
%  
%   IC = setInitialConditions(PDEM, U0) sets the initial conditions for
%   time-dependent problems or the initial-guess for stationary nonlinear 
%   problems. U0, the initial value or guess can be defined as follows:
%    (1) A scalar representing a constant initial value/guess throughout 
%        the domain.
%    (2) A column vector of length PDESystemSize representing a constant 
%        initial value/guess for each equation throughout the domain.
%    (3) A function handle representing a spatially varying initial value/guess 
%        throughout the domain. 
%   IC is an InitialCondition object that associates the initial values with 
%   the geometric domain.
%  
%   IC = setInitialConditions(PDEM, U0, UT0) sets the initial conditions for
%   time-dependent problems. U0 represents the initial value and UT0 the
%   initial time derivative value. U0 and UT0 can be defined as follows:
%    (1) A scalar representing a constant initial value or derivative value 
%        throughout the domain.
%    (2) A column vector of length PDESystemSize representing a constant initial 
%        value or derivative value for each equation throughout the domain.
%    (3) A function handle representing a spatially varying initial value or
%        derivative value throughout the domain. 
%   IC is an InitialCondition object that associates the initial values and
%   derivatives with the geometric domain.
%  
%   IC = setInitialConditions(PDEM, U0, REGIONTYPE, REGIONID) and 
%   IC = setInitialConditions(PDEM, U0, UT0, REGIONTYPE, REGIONID) sets the 
%   initial conditions/guess on a region of the domain defined by REGIONTYPE 
%   and REGIONID. Where REGIONTYPE is 'cell' (3-D only), 'face', 'edge', or 
%   'vertex' and REGIONID is the ID of the geometric region. IC is an 
%   InitialCondition object that associates the initial values and/or 
%   derivatives with the application region.
%
%   IC = setInitialConditions(PDEM, TDRESULTS, Ti) sets the initial conditions 
%   for time-dependent problems or the initial-guess for stationary nonlinear 
%   problems. TDRESULTS is the solution from a previous time-dependent analysis 
%   on the same geometry and mesh. Ti is the time index of TDRESULTS that
%   defines the values to be used for the initial condition. If Ti is not 
%   specified, the last time step is used. IC is an InitialCondition object 
%   that associates the initial values and/or derivatives with the nodes in 
%   the mesh.
%
%   IC = setInitialConditions(PDEM, SRESULTS) sets the initial conditions for
%   time-dependent problems or the initial-guess for stationary nonlinear 
%   problems. SRESULTS is the solution from a previous stationary analysis 
%   on the same geometry and mesh. IC is an InitialCondition object that 
%   associates the initial values and/or derivatives with the nodes in the mesh.
%
%   Precedence rules: 
%   (1) For multiple assignments to the same region, the last assignment wins. 
%   (2) Assignments to lower dimension entities take precedence over assignments 
%       to higher dimensional entities.  For example, if different ICs are 
%       assigned to a face and one of its bounding edges, the edge IC takes 
%       precedence over the face IC.
%   (3) When initial conditions are set using a results object that assignment 
%       supersedes all prior assignments, if any. Subsequent region-based 
%       assignments follow rules (1) and (2).
%  
%   Function handle format: A function handle representation of the ICs is
%   required to accept a single input argument and return a single output 
%   argument. The function is of the form:
%  
%                  Uval = myicfun(locations)
%  
%   The input argument is a struct representing nodal locations within the 
%   assignment region. The fields are: locations.x, locations.y and 
%   for 3-D problems locations.z. The output argument Uval is an MxN matrix
%   where M is PDESystemSize and N is length(locations.x).
%  
%   Example: Set ICs for a 3-D Heat Transfer problem in a block
%      pdem = createpde;
%      importGeometry(pdem,'Block.stl');
%      pdegplot(pdem,'FaceLabels', 'on');
%      % Assign an ambient temperature of 20 to the domain and a local
%      % temperature of 1000 to face 3.
%      ic1 = setInitialConditions(pdem, 20)
%      ic2 = setInitialConditions(pdem, 1000, 'face', 3)
%  
%   See also pde.PDEModel, pde.PDEModel/PDESystemSize

% Copyright 2015-2016 The MathWorks, Inc.


    mcoefdefined= false;
    if ~isempty(self.EquationCoefficients)
        mcoefdefined = self.EquationCoefficients.mDefined();
    end
    narginchk(2,5);             
    nargoutchk(0,1);         
    if isempty(self.Geometry)
         error(message('pde:pdeModel:initialCondNoGeom'));
    end
     
    isnodal = false;
    if isa(varargin{1}, 'pde.PDEResults') 
       if isa(varargin{1}, 'pde.EigenResults') 
          error(message('pde:pdeInitialConditions:invalidResultsObject')); 
       end        
       if numel(varargin) == 1
          argsToPass = {varargin{1},mcoefdefined};        
       elseif numel(varargin) == 2
          argsToPass = {varargin{1},  varargin{2}, mcoefdefined};  
       else
          error(message('pde:pdeModel:setICMultipleVarArgs'));  
       end      
       isnodal = true;
    else 
       argsToPass = ParseGeometricICs(self, varargin{:});
    end
                  
    if isempty(self.InitialConditions)
        iccont = pde.InitialConditionsRecords(self);     
        argsToPass = {iccont, argsToPass{:}};
        if isnodal
            ic = pde.NodalInitialConditions(argsToPass{:}); 
        else
            ic = pde.GeometricInitialConditions(argsToPass{:}); 
        end
        iccont.InitialConditionAssignments = ic; 
        self.InitialConditions = iccont;
        addlistener(iccont,'ObjectBeingDestroyed',@pde.InitialConditionsRecords.preDelete);             
    else      
        argsToPass = {self.InitialConditions, argsToPass{:}};
        if isnodal
            ic = pde.NodalInitialConditions(argsToPass{:}); 
        else
            ic = pde.GeometricInitialConditions(argsToPass{:});
        end
        self.InitialConditions.InitialConditionAssignments(end+1) = ic;         
    end   
    if isnodal
        addlistener(ic,'ObjectBeingDestroyed',@pde.NodalInitialConditions.preDelete);  
    else
        addlistener(ic,'ObjectBeingDestroyed',@pde.GeometricInitialConditions.preDelete);         
    end
      
end

function argsToPass = ParseGeometricICs(self, varargin)
    if ~(isnumeric(varargin{1}) || isa(varargin{1}, 'function_handle'))
             error(message('pde:pdeModel:initialCondBadValue')); % Invalud U0
    end
    if numel(varargin) == 2 || numel(varargin) == 4
        if ~(isnumeric(varargin{2}) || isa(varargin{2}, 'function_handle'))
            error(message('pde:pdeModel:initialCondBadDerivative')); % Invalud UT0 
        end
    end
    inputargs = {};
    if numel(varargin) == 3 
     inputargs = varargin(2:end);
    elseif numel(varargin) == 4 
     inputargs = varargin(3:end);         
    end
         
    parser = inputParser;   
    parser.KeepUnmatched=true;
    parser.addParameter('cell', [], @pde.PDEModel.isValidEntityID);  
    parser.addParameter('face', [], @pde.PDEModel.isValidEntityID); 
    parser.addParameter('edge', [], @pde.PDEModel.isValidEntityID);  
    parser.addParameter('vertex', [], @pde.PDEModel.isValidEntityID);  
    parser.parse(inputargs{:});   
            
    if numel(parser.UsingDefaults) == 4 % Global assignment
      if self.IsTwoD     
           argsToPass = [varargin, {'face', 1:self.Geometry.NumFaces}]; 
      else     
           argsToPass = [varargin, {'cell', 1:self.Geometry.NumCells}]; 
      end 
    else   % Local assignment
      if ~isempty(parser.Results.cell)
        self.isValidEntityID(parser.Results.cell);
        if any(parser.Results.cell > self.Geometry.NumCells)
            error(message('pde:pdeModel:invalidCellIndex'));
        end
      end
      if ~isempty(parser.Results.face)
        self.isValidEntityID(parser.Results.face);
        if any(parser.Results.face > self.Geometry.NumFaces)
            error(message('pde:pdeModel:invalidFaceIndex'));
        end
      end
      if ~isempty(parser.Results.edge)
        self.isValidEntityID(parser.Results.edge);
        if any(parser.Results.edge > self.Geometry.NumEdges)
            error(message('pde:pdeModel:invalidEdgeIndex'));
        end
      end
      if ~isempty(parser.Results.vertex)
        self.isValidEntityID(parser.Results.vertex);
        if any(parser.Results.vertex > self.Geometry.NumVertices)
            error(message('pde:pdeModel:invalidVertexIndex'));
        end
      end  
     argsToPass = varargin;
    end

end

