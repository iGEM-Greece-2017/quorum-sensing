function gm = geometryFromEdges(self, gedges)
% geometryFromEdges       - Create 2D geometry from DECSG or PDEGEOM
%    G = geometryFromEdges(PDEM, DECSG) constructs a geometry object from 
%    DECSG - a decomposed geometry matrix. The geometry created is assigned
%    to the Geometry property of this class and a handle to the geometry is 
%    returned in G.
%
%    G = geometryFromEdges(PDEM, PDEGEOM) constructs a geometry object from 
%    PDEGEOM - a function handle to a geometry file. The geometry created 
%    is assigned to the Geometry property of this class and a handle to the 
%    geometry is returned in G.
%
%    Example: Create geometry from a boundary edge definition in DECSG format
%      First step is to create the geometry description (GD) matrix in CSG  
%      format representing a rectangle with a hole, and then convert to DECSG
%      to get the boundary edges. 
%      In the GD matrix for the rectangle, the first element 3 is the
%      rectangle code, 4 is the number of sides, followed by X and then Y
%      coordinates.
%      R1 = [3,4,-1,1,1,-1,-.4,-.4,.4,.4]';
%      % Now the circle, center (0, 0) radius 0.2
%      C1 = [1,0,0,.2]';
%      % Pad C1 with zeros to enable concatenation with R1
%      C1 = [C1;zeros(length(R1)-length(C1),1)];
%      geom = [R1,C1];
%      % Names for the two geometric objects
%      ns = (char('R1','C1'))';
%      % Set formula
%      sf = 'R1-C1';
%      % Create geometry
%      gd = decsg(geom,sf,ns);
%      % Now create the class geometry representation
%      pdem = createpde();
%      gm = geometryFromEdges(pdem, gd)
%      % Plot the geometry
%      pdegplot(gm, 'EdgeLabels','on')
%      axis equal
%      % Observe that the geometry is appended to the PDEModel
%      pdem
%
% See also  DECSG, PDEGEOM, pde.AnalyticGeometry

% Copyright 2014-2017 The MathWorks, Inc.

        if ~isempty(self.Geometry)                        
          error(message('pde:pdeModel:noAssemSupport')); 
        end
        gm = pde.AnalyticGeometry(gedges);
        self.Geometry = gm;
    end