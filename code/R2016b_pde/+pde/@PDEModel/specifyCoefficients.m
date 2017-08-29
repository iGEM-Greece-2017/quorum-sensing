function coef = specifyCoefficients(self, varargin)  
% specifyCoefficients - Specify all PDE coefficients over a domain or subdomain
%     The PDE toolbox can solve equations of the form:
%  
%              m*d^2u/dt^2 + d*du/dt - div(c*grad(u)) + a*u = f
%  
%     and the corresponding eigenvalue equation of the form:
%  
%                   - div(c*grad(u)) + a*u = lamda*d*u
%     or
%                   - div(c*grad(u)) + a*u = (lamda^2)*m*u
%  
%     The equation to solve is defined in terms of the coefficients m, d, c, 
%     a, f. This method creates an object representing the coefficients 
%     in a domain or subdomain and appends the object to the EquationCoefficients 
%     property.
%  
%     You can define more than one set of equation coefficients. For example, 
%     for a 2-D geometry that has two subdomains and equation coefficients
%     that differ in each subdomain. However, the same coefficient terms must
%     be present in each equation - they just differ in value.
%     Note: Subdomains are not currently supported in 3-D.
%
%     CA = specifyCoefficients(PDEM, 'm', mcoef, 'd', dcoef, 'c', ccoef,...
%                              'a', acoef, 'f', fcoef)
%     Specifies all of the coefficients of the equations throughout the analysis
%     domain using name-value pairs to define them. If a coefficient is not 
%     defined in the PDE, specify 0. CA is a CoefficientAssignment object that 
%     associates the equation with the geometric domain.
%     The coefficient names are as follows:
%  
%     Name   Coefficient Of                  Typical Physical Effect
%     --------------------------------------------------------------
%     'm'    Second-order derivative         Mass
%     'd'    First-order derivative          Damping or mass
%     'c'    Grad term                       Diffusion
%     'a'    Dependent variable u            Absorption
%     'f'    Right-hand side                 Source term
%
%     The corresponding Value of the coefficients can be defined in one of two ways: 
%        (1) - Numeric format as a scalar, vector or matrix
%        (2) - MATLAB function format
%
%     Note: If the PDE to solve has both 'm' and 'd' coefficients that are
%           non-zero, then 'd' must be a matrix that has the same size as
%           the discretized mass matrix (not to be confused with 'm').
%     
%     CA = specifyCoefficients(..., 'face', FACEIDS)
%     Specifies all of the coefficients of the equations defined within specific
%     faces of a 2-D Geometric Model. FACEIDS is a vector of face IDs, 
%     where each 0 < FACEIDS(j) < NumFaces in the geometry. CA is a 
%     CoefficientAssignment object that associates the equation with the faces 
%     that define a subdomain.
%
%     Example 1: Assign the same equation coefficients throughout the 
%        % domain of a 3-D model.
%        pdem = createpde();
%        importGeometry(pdem, 'Block.stl');
%        pdegplot(pdem);
%        ca = specifyCoefficients(pdem,'m',0,'d',1,'c',1,'a',0,'f',1)
%
%     Example 2: Assign different equation coefficients in each of the 
%        % subdomains of a 2-D model.
%        pdem = createpde();
%        g = @lshapeg;
%        geometryFromEdges(pdem,g);
%        pdegplot(pdem, 'subdomainLabels','on')
%        ca1 = specifyCoefficients(pdem,'m',0,'d',1,'c',1,'a',0,'f',1,'face',1)
%        ca2 = specifyCoefficients(pdem,'m',0,'d',1,'c',1,'a',0,'f',2,'face',2)
%        ca3 = specifyCoefficients(pdem,'m',0,'d',1,'c',1,'a',0,'f',3,'face',3)
%        
% 
% See also pde.PDEModel, assembleFEMatrices, NumericCoefficientFormat, FunctionCoefficientFormat

% Copyright 2015 The MathWorks, Inc.

    narginchk(11,13); 
    nargoutchk(0,1); 
    if nargin == 12
         error(message('pde:pdeModel:invalidNumInputArgs'));  
    end
     
    if isempty(self.Geometry)
         error(message('pde:pdeModel:coefAssignNoGeom'));
    end
    
    argsToPass = [varargin, {'SystemSize', self.PDESystemSize}];
    parser = inputParser;
    parser.PartialMatching=false; % Clash between 'face' and 'f'     
    parser.addParameter('face', [], @pde.PDEModel.isValidEntityID); 
    parser.addParameter('m', []);                 
    parser.addParameter('d', []);
    parser.addParameter('c', []);
    parser.addParameter('a', []);    
    parser.addParameter('f', []); 
    
    parser.parse(varargin{:});
    
    if ~isempty(parser.UsingDefaults) && ~(strcmp(parser.UsingDefaults{1},'face'))
        error(message('pde:pdeModel:missingSpecifyCoefficients'));
    end
    
    if ~isempty(parser.UsingDefaults) && strcmp(parser.UsingDefaults{1},'face')
      if self.IsTwoD
          argsToPass = [{'face', 1:self.Geometry.NumFaces}, varargin, {'SystemSize', self.PDESystemSize}]; 
        else
          argsToPass = [{'cell', 1:self.Geometry.NumCells}, varargin, {'SystemSize', self.PDESystemSize}]; 
      end 
    else
      if ~self.IsTwoD  
          error(message('pde:pdeModel:noCoefOnBoundary')); 
      end
      self.isValidEntityID(parser.Results.face);
      if any(parser.Results.face > self.Geometry.NumFaces)
          error(message('pde:pdeModel:invalidFaceIndex'));
      end
    end
      
    if isempty(self.EquationCoefficients)
        coefcont =  pde.CoefficientAssignmentRecords(self);             
        coef = pde.CoefficientAssignment(coefcont,argsToPass{:});      
        coefcont.CoefficientAssignments = coef; 
        self.EquationCoefficients = coefcont;
        addlistener(coefcont,'ObjectBeingDestroyed',@pde.CoefficientAssignmentRecords.preDelete);       
    else
        coef = pde.CoefficientAssignment(self.EquationCoefficients,argsToPass{:});
        self.EquationCoefficients.CoefficientAssignments(end+1) = coef;         
    end   
    addlistener(coef,'ObjectBeingDestroyed',@pde.CoefficientAssignment.preDelete);         
      
end