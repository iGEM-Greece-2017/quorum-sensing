function F = facetAnalyticGeometry(~, ag)
%facetAnalyticGeometry - Facet 2D geometry in DECSG or PDEGEOM format.
%
% Internal use only

%   Copyright 2017 The MathWorks, Inc.
    g = ag.geom;
    nbs=pdeigeom(g); % Number of boundary segments
    d=pdeigeom(g,1:nbs); % Row 1,2 are start and end params, 3&4 L/R faces
    NSEGPERCURVE=12;   
    % Xg, Yg, Cg arrays for geometry
    Xg = zeros(0,1);
    Yg = zeros(0,1);
    Cg = zeros(0,2);
    CidToEid = zeros(0,1);
    LRFaces = d(3:4,:)';
    for i = 1:nbs         
       % sample the curve 
       params = linspace(d(1, i), d(2, i), NSEGPERCURVE+1); 
       params(1) = d(1, i);
       params(end) = d(2, i);      
       [x, y] = pdeigeom(g,i, params); 
       startn = numel(Xg)+1;
       endn = numel(Xg) + numel(x);
       c = [(startn:(endn-1))' (startn+1:endn)'];
       nump = numel(x);
       Xg(end+(1:nump)) = x;
       Yg(end+(1:nump)) = y; 
       Cg(end+(1:nump-1),:) = c;
       CidToEid(end+(1:nump-1),:) = i;
    end


    % Merge duplicate points and rewire constraint IDs
    comptol = max(max(abs(Xg)), max(abs(Yg)));
    comptol = max(comptol,1)*eps*1000;     
    [~, I, IC] = uniquetol([Xg, Yg], comptol, 'ByRows',true);
    ca = Cg(:,1);
    cb = Cg(:,2);
    Cg = [I(IC(ca)) I(IC(cb))];
    I = sort(I);
    ne = numel(I);
    mymap = containers.Map(I, 1:ne);
    ne = numel(Cg);
    for i = 1:ne
       Cg(i) = mymap(Cg(i));
    end
    P = [Xg(I) Yg(I)];


    numCpre = size(Cg,1);
    warnState(1) = warning('off','MATLAB:delaunayTriangulation:ConsSplitPtWarnId');
    warnState(2) = warning('off', 'MATLAB:delaunayTriangulation:ConsConsSplitWarnId');
    warnState(3) = warning('off', 'MATLAB:delaunayTriangulation:DupConsWarnId');
    warnState(4) = warning('off', 'MATLAB:delaunayTriangulation:LoopConsWarnId');
    warnState(5) = warning('off','MATLAB:delaunayTriangulation:DupPtsWarnId');
    warnState(6) = warning('off','MATLAB:delaunayTriangulation:DupPtsConsUpdatedWarnId');
    dt = delaunayTriangulation(P, Cg);
    warning(warnState);
    numCpost = size(dt.Constraints,1);
    if numCpre ~= numCpost
       error(message('pde:pdeModel:BadGeomIntersectingEdges'));      
    end    
    % Check for a spur edges
    % Each vertex should have a valency of two or more
    CC = dt.Constraints;
    [valency,~]=hist(CC(:),unique(CC(:)));
    if any(valency == 1)
       error(message('pde:pdeModel:BadGeomSpurEdges'));   
    end    
    F.Points = [P zeros(size(P,1),1)];
    F.Triangles = dt.ConnectivityList;
    F.Constraints = Cg;
    vxy = ag.vertexLabelLocations();
    dt = delaunayTriangulation(P);
    F.VxIdToPointId = dt.nearestNeighbor(vxy{1}, vxy{2});
    F.ConsIdToEdgeId = CidToEid;
    F.LRFaces = LRFaces;
end