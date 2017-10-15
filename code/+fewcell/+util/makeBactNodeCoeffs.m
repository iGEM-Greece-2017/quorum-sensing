function makeBactNodeCoeffs(mesh,nBact,geom)
% Find the mesh nodes that correspond to each bacterium in the geometry and calculate the 
% coefficient that translates the singlecell dAHL into a change for the corresponding mesh nodes

global bactNodes;
global bactNodeIdx;
global bactNodesEqulength;
printPerBact= nBact < 9;
faceIdx= 4;
if strcmp(mesh.GeometricOrder, 'quadratic'), faceIdx= 7; end
[p,~,t]= meshToPet(mesh);
%% Order bacterial subdomains
% Scan the geometry in a radial-out direction (all layers per ring) and assign each bacterium to its
% corresponding subdomain number
bactSubdomain= zeros(nBact,1);
if nargin==2
  bactSubdomain= 2:nBact+1;
else
  % Upper-left-corner node coordinates
  ybact= -(0:geom.nLayers-1)'*(geom.bactSize(2)+geom.layerSeparation);
  xbact= geom.bactCenter0(1)-geom.bactSize(1)/2 + (0:geom.nRings-1)'*(geom.bactSize(1)+geom.ringDist);
  [X,Y]= meshgrid(xbact,ybact);
  bactCoord= [X(:)';Y(:)'];
  % Find X,Y idx in p
  for b= 1:nBact
    idx= find(sum(abs(bactCoord(:,b)-p)<1e-7)==2);
    bactSubdomain(b)= setdiff(unique( t(faceIdx,any(t(1:end-1,:)==idx)) ),1);
  end
end
%% 
%{
% Find the mesh nodes that correspond to each bacterium in the geometry
bactNodes= cell(nBact,1);
% Find bacterial nodes
for b= 1:nBact
  bactNodes{b}= unique( t(1:end-1,t(faceIdx,:)==b+1) );         % all nodes for elements in bacterium
  if isempty(bactNodes{b}), error(['No mesh nodes on bacterium #',num2str(b)]); end
  if printPerBact, fprintf('Nodes in bact#%d: %d\n', b,length(bactNodes{b})); end
end

if ~printPerBact
  bactNodeN= cellfun(@length, bactNodes);
  fprintf('[makeBactNodeCoeffs]: Too many bacteria to list. Printing max/min node counts\n');
  fprintf('Min bact nodes: %d\nMax bact nodes: %d\n', min(bactNodeN), max(bactNodeN));
end
%}

bactNodes= sparse(size(p,2), nBact, false);
occupied= sparse(size(p,2),1, false);
for b=1:nBact
  candidates= sparse(unique( t(1:end-1,t(faceIdx,:)==bactSubdomain(b)) ), 1, true);
  bactNodes(candidates & ~occupied(1:length(candidates)), b)= true;
  occupied(bactNodes(:,b))= true;
  if ~nnz(bactNodes(:,b)), error(['No mesh nodes on bacterium #',num2str(b)]); end
  if printPerBact, fprintf('Nodes in bact#%d: %d\n', b,nnz(bactNodes(:,b))); end
end
bactNodesEqulength= length(unique(sum(bactNodes,1))) == 1;

bactNodeIdx= [];
if bactNodesEqulength
  fprintf('[makeBactNodeCoeffs]: All bacteria have the same number of nodes. Optimization enabled\n');
  [bactNodeIdx,~]= find(bactNodes);
  bactNodeIdx= reshape(bactNodeIdx, [], nBact);
end

if ~printPerBact
  bactNodeN= full(sum(bactNodes,1));
  fprintf('[makeBactNodeCoeffs]: Too many bacteria to list. Printing max/min node counts\n');
  fprintf('Min bact nodes: %d\nMax bact nodes: %d\n', min(bactNodeN), max(bactNodeN));
end
