classdef StationaryResults < pde.PDEResults
    % pde.StationaryResults PDE solution and its derived quantities
    %   A StationaryResults object provides a convenient representation of
    %   results data for stationary PDE analysis, together with
    %   interpolation and gradient evaluation functions. Create a
    %   StationaryResults object using the createPDEResults function.
    %
    % StationaryResults methods:
    %   interpolateSolution  - Interpolate solution at specified spatial
    %                          locations
    %   evaluateGradient     - Evaluate gradients of solution at specified
    %                          spatial locations
    %   evaluateCGradient    - Evaluate tensor product of c-coefficient and
    %                          gradients of the PDE solution
    %
    % StationaryResults properties:
    %   NodalSolution        - Solution to PDE at nodal locations
    %   XGradients           - Spatial gradient of solution along x-direction.
    %   YGradients           - Spatial gradient of solution along y-direction.
    %   ZGradients           - Spatial gradient of solution along z-direction.
    %   Mesh                 - Discretization of the domain
    %
    %    See also createPDEResults
    
    % Copyright 2015-2016 The MathWorks, Inc.
    
    
    properties(SetAccess = protected)
        % NodalSolution - Solution to PDE at nodal locations
        % Solution array which can be conveniently indexed into to extract
        % a sub-array of interest. The shape of NodalSolution depends on
        % the type of PDE and solver settings. It will be a:
        %
        %   column vector - for a single PDE with no time dependency
        %   matrix        - for a single hyperbolic or parabolic problem,
        %                   or a system of elliptic problems, or a single
        %                   eigenvalue problem
        %   3-D array     - for a system of hyperbolic, parabolic, or
        %                   eigenvalue problems
        %
        % The first array dimension of NodalSolution represents node index.
        % The second array dimension represents the time-step or
        % eigenvector index for a single PDE, or the equation index for a
        % system of PDEs. The third array dimension represents the
        % time-step index for a system of time-dependent PDEs, or the
        % eigenvect index for an eigenvalue problem involving a system of
        % PDEs.
        NodalSolution;
        
        % XGradients - Spatial gradient of solution along x-direction. The
        % shape of the XGradients array is identical to NodalSolution.
        XGradients;
        % YGradients - Spatial gradient of solution along y-direction. The
        % shape of the YGradients array is identical to NodalSolution.
        YGradients;
        % ZGradients - Spatial gradient of solution along z-direction. The
        % shape of the ZGradients array is identical to NodalSolution.
        ZGradients;
    end
    
    methods
      function this= hackitConstructor(this,args)
        u= args{3}{1}; t= args{3}{2};
        this.IsTimeEig= 0;
        this.NodalSolution= u.NodalSolution(:,t);
        this.XGradients= u.XGradients(:,t); this.YGradients= u.YGradients(:,t);
        
        [p,~,t] = meshToPet(u.Mesh);
        this.Interpolant= pdeInterpolant(p,t,this.NodalSolution);
        this.InterpolantdUdx= pdeInterpolant(p,t,this.XGradients);
        this.InterpolantdUdy= pdeInterpolant(p,t,this.YGradients);
      end
      
        function obj = StationaryResults(varargin)
            obj@pde.PDEResults(varargin{:});
            if nargin == 0
                return
            end
            if ischar(varargin{nargin})
              if strcmp(varargin{nargin}, 'hackit, please!')
                obj= obj.hackitConstructor(varargin);
                return;
              end
            end
            narginchk(2,2);
            pdem = varargin{1};
            u = varargin{2};
            if (obj.IsTimeEig)
                error(message('pde:PDEResults:notStationaryResults'));
            end
            ureshaped = pde.PDEResults.reshapePDESolution(u, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            obj.NodalSolution = ureshaped;
            if ~isempty(pdem.EquationCoefficients)
                obj.Coefficients = pdem.EquationCoefficients.packCoefficients;
            else
                obj.Coefficients = [];
            end
            
            obj.XGradients = [];
            obj.YGradients = [];
            obj.ZGradients = [];
            obj.InterpolantdUdx = [];
            obj.InterpolantdUdy = [];
            obj.InterpolantdUdz = [];
            if isempty(obj.NodalSolution)
                return;
            end
            
            obj.Interpolant   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                ureshaped, pdem.PDESystemSize, obj.NumTimeEig);
            
            % Calculate gradients, assign nodal gradients as properties and
            % construct interpolants for gradients
            [ux,uy,uz] = nodalGradients(obj.Interpolant);
            obj.XGradients = pde.PDEResults.reshapePDESolution(ux, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            obj.YGradients = pde.PDEResults.reshapePDESolution(uy, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            obj.ZGradients = pde.PDEResults.reshapePDESolution(uz, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            
            obj.InterpolantdUdx   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.XGradients, pdem.PDESystemSize, obj.NumTimeEig);
            
            obj.InterpolantdUdy   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.YGradients, pdem.PDESystemSize, obj.NumTimeEig);
            
            if ~obj.IsTwoD
                obj.InterpolantdUdz   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                    obj.ZGradients, pdem.PDESystemSize, obj.NumTimeEig);
            end
            
            
        end
        
        % Methods declaration
        function uintrp = interpolateSolution(obj,varargin)
          uintrp = pde.PDEResults.interpolateSolutionInternal(obj,varargin{:});
        end
        function [dudxInt, dudyInt, dudzInt] = evaluateGradient(obj,varargin)
          if (isempty(obj.InterpolantdUdx))
              error(message('pde:PDEResults:gradientNotSupported15bObj'));
          end
          [dudxIntrp, dudyIntrp, dudzIntrp] = pde.PDEResults.evalGradInternal(obj,varargin{:}) ;
        end
        function [cdudx, cdudy, cdudz] = evaluateCGradient(obj,varargin)
          cCoeffAtNodes = pde.PDEResults.evalCCoeffAtNodes(obj);
          dudx = obj.XGradients';
          dudy = obj.YGradients';
          if obj.IsTwoD
              % 2-D implementation code
              [cdudx, cdudy] = pde.PDEResults.cgradImpl2D(obj,cCoeffAtNodes,dudx,dudy);
              cdudz = [];
          else
              dudz = obj.ZGradients';
              % 3-D implementation code
              [cdudx, cdudy, cdudz] = pde.PDEResults.cgradImpl3D(obj,cCoeffAtNodes,dudx,dudy,dudz);
          end
          % Perform interpolation if user specifies required additional arguments.
          if (nargin > 1)
              [cdudx, cdudy, cdudz] = pde.PDEResults.evalCGradInternal(obj,cdudx,cdudy,cdudz,varargin{:});
          end
        end
    end
    
    
    methods (Hidden = true, Static = true, Access = protected)
        function obj = loadobj(obj)
            % function called during loading an object of this type from
            % a MAT-file, interpolant objects need to be constructed
            obj.Interpolant = pde.PDEResults.constructInterpolat(obj.Mesh, ...
                obj.NodalSolution, obj.PDESystemSize, obj.NumTimeEig);
            
            obj.InterpolantdUdx   = pde.PDEResults.constructInterpolat(obj.Mesh,...
                obj.XGradients, obj.PDESystemSize, obj.NumTimeEig);
            
            obj.InterpolantdUdy   = pde.PDEResults.constructInterpolat(obj.Mesh,...
                obj.YGradients, obj.PDESystemSize, obj.NumTimeEig);
            
            if ~obj.IsTwoD
                obj.InterpolantdUdz   = pde.PDEResults.constructInterpolat(obj.Mesh,...
                    obj.ZGradients, obj.PDESystemSize, obj.NumTimeEig);
            end
        end
        
    end
    
    
    
end
