function makeBactNodeCoeffs(mesh,bactSize,bactRingDensity)
% Find the mesh nodes that correspond to each bacterium in the geometry and calculate the 
% coefficient that translates the singlecell dAHL into a change for the corresponding mesh nodes

global bactNodes;
global bactNodeCoeffs;
[p,~,t]= meshToPet(mesh);

% Find the mesh nodes that correspond to each bacterium in the geometry
nBact= length(unique(t(4,:)))-1;
bactNodes= cell(nBact,1);
% Find bacterial nodes
domainNodes= unique( t(1:3,t(4,:)==1) );
for b= 1:nBact
  bactNodes{b}= unique( t(1:3,t(4,:)==b+1) );         % all nodes for elements in bacterium
  internalNodes= setdiff(bactNodes{b}, domainNodes);
  % strip boundary nodes
  if ~isempty(internalNodes), bactNodes{b}= internalNodes; end
  if isempty(bactNodes{b}), error(['No mesh nodes inside bacterium #',num2str(b)]); end
  fprintf('Nodes in bact#%d: %d\n', b,length(bactNodes{b}));
end

% Calculate dAHL mesh node coefficients
% Each node is part of many triangles. Each triangle has an area relative to the bacterium
% Each node's relative area is the sum of the relative areas of each triangle it is part of
% divided by the number of nodes for each triangle (3)
bactNodeCoeffs= cell(nBact,1);
assert(strcmp(mesh.GeometricOrder,'linear'));
nodePerElt= 3;    % because the mesh is linear triangular
for b= 1:nBact
  selElt= t(4,:)==b+1;
  nodesForeachElt= t(1:3,selElt);
  eltRelArea= pdetrg(p,nodesForeachElt)./prod(bactSize);
  nodes= bactNodes{b};
  nodeRelArea= zeros(length(nodes),1);
  for n= 1:length(nodes)
    nodePartofElt= any(nodesForeachElt == nodes(n));
    nodeRelArea(n)= sum(eltRelArea(nodePartofElt))./nodePerElt;
  end
  bactNodeCoeffs{b}= bactRingDensity(b).*nodeRelArea;
end

%ahlProd= dy(6).*bactRingDensity(b).*([1;2;1;2;3;3]./3*0.25);  %<dydt>*<elt/node>/<nodes/elt>*<eltRelArea>
