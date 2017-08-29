function bc = applyBoundaryCondition(self,varargin) 
% applyBoundaryCondition - Apply a boundary condition (BC) to the geometry.  
%      Creates a BC object that applies a BC to the geometry stored in the 
%      Geometry property. The BC object is appended to the BoundaryConditionsRecords 
%      and a handle to the BC is returned. 
% 
%      The PDE toolbox supports Dirichlet BC of the form:
%    
%                h*u = r
%      
%      and Generalized Neumann BC of the form:
% 							
%                n*(c*grad(u)) + q*u = g
%   
%     BC = applyBoundaryCondition(PDEM, BCTYPE,...) sets the boundary conditions 
%     of BCTYPE on a boundary of the domain. BCTYPE is either 'dirichlet',
%     'neumann', or 'mixed'. BCTYPE is an option characters passed in as a positional 
%     argument. BC is a BoundaryCondition object that associates the boundary 
%     values with the application region.
% 
%     BC = applyBoundaryCondition(..., REGIONTYPE, REGIONID,...) sets the boundary
%     conditions on a region of the domain defined by REGIONTYPE and REGIONID. 
%     Where name/value pair REGIONTYPE is 'Face','Edge', 'Vertex' and REGIONID 
%     is an integer specifying the ID of the geometric entity. 
% 
%     BC = applyBoundaryCondition(...,'Name',Value)
%     sets boundary conditions using the following name-value pairs.
% 
%     For Dirichlet boundary conditions, specify either both arguments 'r' 
%     and 'h', or the argument'u'. When specifying 'u', you can also use 'EquationIndex'.
%               'r'             r-coefficient for a Dirichlet boundary condition. 
%                               Vector of N numerical values or a function handle.
%               'h'             h-coefficient for a Dirichlet boundary condition. 
%                               N X N matrix of numerical values or a function 
%                               handle. Default is eye(N).
%               'u'             Vector of N numerical values or a function 
%                               handle. A more compact way of defining the 
%                               most common types of Dirichlet boundary conditions.
%               'EquationIndex' Vector with N integers in the range 1:N defining 
%                               an equation index for each entry in the 'u' 
%                               vector.
% 
%     For Neumann boundary conditions:
%               'g'             g-coefficient for a Neumann boundary condition.
%                               Vector of N numerical values or a function 
%                               handle. Default is zeros(N,1).
%               'q'             q-coefficient for a Neumann boundary condition. 
%                               N X N matrix of numerical values or a function 
%                               handle. Default is zeros(N,N).
% 
%     For mixed boundary conditions, you can use name-value pairs from both 
%     Dirichlet and Neumann boundary conditions as needed. The default boundary 
%     condition for a region of boundary is zero Neumann, with g = 0.
% 
%               'Vectorized'    If set to 'on', indicates that the custom
%                               function can accept input arguments at multiple
%                               locations and returns boundary condition
%                               matrices defined at these multiple locations.
%                               The default is 'off'. This setting is useful when one of the
%                               other parameters are defined as a function handle.
% 
%     Example 1: Set the value of the unknown, u, to zero on a single edge
%            pdem = createpde(1);
%            R1 = [3,4,-1,1,1,-1,-.4,-.4,.4,.4]';
%            g = decsg(R1);
%            pg = geometryFromEdges(pdem,g);
%            NumEdges = pg.NumEdges;
%            bc1 = applyBoundaryCondition(pdem,'dirichlet','edge',1,'u',0);
% 
%     Example 2: Set the value of the unknown, u, to two on three edges
%            pdem = createpde(1);
%            R1 = [3,4,-1,1,1,-1,-.4,-.4,.4,.4]';
%            g = decsg(R1);
%            pg = geometryFromEdges(pdem,g);
%            NumEdges = pg.NumEdges;
%            bc1 = applyBoundaryCondition(pdem,'dirichlet','edge',(1:3),'u',2);
% 
%     Example 3: Set the value of u(2) to three on edge four
%            numberOfPDEs = 2;
%            pdem = createpde(numberOfPDEs);
%            R1 = [3,4,-1,1,1,-1,-.4,-.4,.4,.4]';
%            g = decsg(R1);
%            pg = geometryFromEdges(pdem,g);
%            NumEdges = pg.NumEdges;
%            bc1 = applyBoundaryCondition(pdem,'mixed','edge',4,'u',3,'EquationIndex',2);
% 
%     Example 4: Same boundary condition as example 3 using 'r' and 'h'
%            numberOfPDEs = 2;
%            pdem = createpde(numberOfPDEs);
%            R1 = [3,4,-1,1,1,-1,-.4,-.4,.4,.4]';
%            g = decsg(R1);
%            pg = geometryFromEdges(pdem,g);
%            NumEdges = pg.NumEdges;
%            h = [0 0;0 1];
%            r = [0 3];
%            bc1 = applyBoundaryCondition(pdem,'mixed','edge',4,'r',r,'h',h);
% 
%     Example 5: Set Neumann g vector 
%            numberOfPDEs = 2;
%            pdem = createpde(numberOfPDEs);
%            R1 = [3,4,-1,1,1,-1,-.4,-.4,.4,.4]';
%            g = decsg(R1);
%            pg = geometryFromEdges(pdem,g);
%            NumEdges = pg.NumEdges;
%            bc1 = applyBoundaryCondition(pdem,'neumann','edge',4,'g',[0 .123]);
% 
%     Example 6: Set Mixed BCs in a block for N = 3
%            pdem = createpde(3);
%            importGeometry(pdem,'Block.stl');
%            % Apply to face 3 generalized Neumann BC for equation 1 and 
%            % Dirichlet BC or equation 2 and 3.
%            h = [0 0 0;0 1 0;0 0 1];
%            r = [0 3 3];
%            q = [1 0 0;0 0 0;0 0 0];
%            g = [3 0 0];
%            bc = applyBoundaryCondition(pdem,'mixed','face',3,'h',h,'r',r,'g',g,'q',q);
% 
%     Example 7: Same boundary condition as example 6 using 'u' and 'EquationIndex'
%            pdem = createpde(3);
%            importGeometry(pdem,'Block.stl');
%            u = [3 3];
%            q = [1 0 0;0 0 0;0 0 0];
%            g = [3 0 0];
%            bc = applyBoundaryCondition(pdem,'mixed','face',3,...
%             'u',u,'EquationIndex',[2 3],'g',g,'q',q);
% 
%        See also pde.BoundaryCondition, pde.PDEModel/geometryFromEdges, DECSG
% 
%        Copyright 2014-2016 The MathWorks, Inc.

    narginchk(3,12); 
    nargoutchk(0,1);
    
    argsToPass = [varargin,{'SystemSize', self.PDESystemSize}];
    
    nc = numel(varargin{1});
    if ~(strncmp(varargin{1},'dirichlet',nc) || strncmp(varargin{1},'neumann',nc) || strncmp(varargin{1},'mixed',nc)...
            || strncmpi(varargin{1},'face',nc) || strncmpi(varargin{1},'edge',nc) || strncmpi(varargin{1},'r',nc)...
            || strncmpi(varargin{1},'h',nc) || strncmpi(varargin{1},'u',nc) || strncmpi(varargin{1},'g',nc)...
            || strncmpi(varargin{1},'q',nc) || strncmpi(varargin{1},'EquationIndex',nc) || strncmpi(varargin{1},'Vectorized',nc))
        % not a valid first argument
        error(message('pde:pdeBoundaryConditions:invalidBCType'));
    end
    
    if ~(strncmp(varargin{1},'dirichlet',nc) || strncmp(varargin{1},'neumann',nc) || strncmp(varargin{1},'mixed',nc))
         % this must be an old syntax
        argsToPass = [{'unknown'},argsToPass];
    end
    
    argsToTest = argsToPass(2:end-2);
   
    if ~isempty(self.Geometry) && numel(varargin) > 1  
        parser = inputParser;                   
        parser.KeepUnmatched=true;  
        parser.addParameter('Face', [], @isreal);
        parser.addParameter('Edge', [], @isreal);        
        parser.parse(argsToTest{:}); 
        if ~isempty(parser.Results.Face) 
            validateattributes(parser.Results.Face,{'numeric'},{'integer','positive', 'nonzero','real', 'nonsparse'});           
            if any(parser.Results.Face > self.Geometry.NumFaces)
              error(message('pde:pdeModel:invalidFaceIndex'));
            end
        elseif ~isempty(parser.Results.Edge)           
            validateattributes(parser.Results.Edge,{'numeric'},{'integer','positive','nonzero','real','nonsparse'});           
            if any(parser.Results.Edge > self.Geometry.NumEdges)
              error(message('pde:pdeModel:invalidEdgeIndex'));
            end
        end
    end
    
    if isempty(self.BoundaryConditions)
        BCcont = pde.BoundaryConditionRecords(self);             
        bc = pde.BoundaryCondition(BCcont,argsToPass{:});      
        BCcont.BoundaryConditionAssignments = bc; 
        self.BoundaryConditions = BCcont;
    else
        bc = pde.BoundaryCondition(self.BoundaryConditions,argsToPass{:});
        self.BoundaryConditions.BoundaryConditionAssignments(end+1) = bc;         
    end      
    
end