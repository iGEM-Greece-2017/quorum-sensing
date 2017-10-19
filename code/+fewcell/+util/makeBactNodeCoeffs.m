function bactSubdomain= makeBactNodeCoeffs(mesh,nBact,geom)
% Find the mesh nodes that correspond to each bacterium in the geometry and calculate the 
% coefficient that translates the singlecell dAHL into a change for the corresponding mesh nodes
  global bactNodes;
  global bactNodeIdx;
  global bactNodesEqulength;
  faceIdx= 4;
  if strcmp(mesh.GeometricOrder, 'quadratic'), faceIdx= 7; end
  [p,~,t]= meshToPet(mesh);
  %% Order bacterial subdomains
  % Scan the geometry in a radial-out direction (all layers per ring first) and assign each bacterium to its
  % corresponding subdomain number
  bactSubdomain= zeros(nBact,1);
  if nargin==2
    bactSubdomain= 2:nBact+1;
  else
    % Upper-left-corner node coordinates
    ybact= -(0:geom.nLayers-1)'*(geom.bactSize(2)+geom.layerSeparation*geom.bactSize(2));
    xbact= geom.bactCenter0(1)-geom.bactSize(1)/2 + (0:geom.nRings-1)'*(geom.bactSize(1)+geom.ringDist*geom.bactSize(1));
    [X,Y]= meshgrid(xbact,ybact);
    bactCoord= [X(:)';Y(:)'];
    % Find X,Y idx in p
    for b= 1:nBact
      idx= find(sum(abs(bactCoord(:,b)-p)<1e-7)==2);
      bactSubdomain(b)= setdiff(unique( t(faceIdx,any(t(1:end-1,:)==idx)) ),1);
    end
  end

  %% Find bacterial nodes
  bactNodes= sparse(size(p,2), nBact, false);
  occupied= sparse(size(p,2),1, false);
  for b=1:nBact
    candidates= sparse(unique( t(1:end-1,t(faceIdx,:)==bactSubdomain(b)) ), 1, true);
    bactNodes(candidates & ~occupied(1:length(candidates)), b)= true;
    occupied(bactNodes(:,b))= true;
    if ~nnz(bactNodes(:,b)), error(['No mesh nodes on bacterium #',num2str(b)]); end
  end
  bactNodesEqulength= length(unique(sum(bactNodes,1))) == 1;

  bactNodeIdx= [];
  if bactNodesEqulength
    [bactNodeIdx,~]= find(bactNodes);
    bactNodeIdx= reshape(bactNodeIdx, [], nBact);
  else
    fprintf('[makeBactNodeCoeffs]: All bacteria don''t have the same number of nodes. Optimization disabled\n');
  end

  bactNodeN= full(sum(bactNodes,1));
  if min(bactNodeN) ~= max(bactNodeN)
    fprintf('Min/Max bact nodes: %d - %d\n', min(bactNodeN), max(bactNodeN));
  end
end
