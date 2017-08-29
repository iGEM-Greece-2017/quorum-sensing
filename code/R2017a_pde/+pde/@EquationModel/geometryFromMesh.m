function [g, m] = geometryFromMesh(self, nodes, elements)   
% geometryFromMesh - Creates a 3D geometric model from a triangulated mesh
%    geometryFromMesh Creates a 3D geometric model from the boundary of a 
%    tetrahedral mesh or general surface triangulation. For mesh input, the 
%    mesh is retained by default and can be deleted/regenerated as necessary.
%    
%    geometryFromMesh(PDEM, Nodes, Elements) constructs a geometric model from 
%    a mesh defined by Nodes and Elements. The variable Nodes is a matrix of
%    size 3-by-nnodes, where nnodes is the number of nodes in the mesh. The 
%    variable Elements is a matrix of size mnpe-by-nelem where nelem is the 
%    number of elements in the mesh and mnpe is the number of nodes per 
%    element; 3 for a surface triangulation, 4 or 10 for a volume mesh. The 
%    geometry created is assigned to the Geometry property of this class. 
%    If the input is a volume mesh, it is likewise assigned to the Mesh 
%    property.
%  
%    [G, MSH] = geometryFromMesh(__) returns handles to the geometry G and 
%    mesh MSH. If the input mesh is not a volume mesh, MSH is empty, the mesh 
%    can be generated using the generateMesh method.
%  
%    Example 1: Construct a geometric block from the convex hull of a
%               meshgrid of points.
%      [x,y,z] = meshgrid(-2:4:2);
%      x = x(:); y = y(:); z = z(:);
%      K = convhull(x,y,z);
%      pdem = createpde();
%      geometryFromMesh(pdem, [x'; y'; z'], K')
%      pdegplot(pdem);
%  
%    Example 2: Create block with a cylindrical hole using an alphaShape.
%              Then create a geometry from the boundary triangulation 
%              of the alphaShape.
%      Steps: Create a square grid of points and a radial grid of points
%              contained within it. Boolean these points using an alpha
%              shape, then sweep the result to 3D. 
%      [xg, yg] = meshgrid(-3:0.25:3);
%      xg = xg(:);
%      yg = yg(:);
%      t = (pi/24:pi/24:2*pi)';
%      x = cos(t);
%      y = sin(t);
%      circShp = alphaShape(x,y,2);
%      in = inShape(circShp, xg,yg);
%      xg = [xg(~in); cos(t)];
%      yg = [yg(~in); sin(t)];
%      zg = ones(numel(xg),1);
%      % Sweep the points in Z
%      xg = repmat(xg,5,1);
%      yg = repmat(yg,5,1);
%      zg = zg*(0:.25:1);
%      zg = zg(:);
%      shp = alphaShape(xg,yg,zg, 0.5);
%      [tri, pts] = boundaryFacets(shp);
%      pdem = createpde();
%      geometryFromMesh(pdem,pts', tri');
%      pdegplot(pdem);   
%      generateMesh(pdem);
%      figure
%      pdemesh(pdem)
% 
%    Example 3: Create a geometric model from a volume tetrahedral mesh
%      load tetmesh 
%      % This loads variables X and tet into the workspace
%      % X represents the nodes and tet represents the elements
%      pdem = createpde();
%      geometryFromMesh(pdem, X', tet')
%      pdegplot(pdem);  
%     
%   See also pde.FEMesh, pde.DiscreteGeometry, alphaShape, convhull, delaunay.    

% Copyright 2015-2016 The MathWorks, Inc.


    if ~isempty(self.Geometry)                        
           error(message('pde:pdeModel:noAssemSupport')); 
    end
    try       
        [g, m] = geomFromMeshInternal(nodes, elements);        
    catch ex
        throwAsCaller(ex);
    end 
    self.Geometry = g;
    self.Mesh = m;
    self.IsTwoD = false; 
    if nargout == 2 && isempty(m)
        warning(message('pde:pdeModel:invalidDomainMesh')); 
    end
end
