function makeBactNodeCoeffs(mesh,bactSize,bactRingDensity)
% Find the mesh nodes that correspond to each bacterium in the geometry and calculate the 
% coefficient that translates the singlecell dAHL into a change for the corresponding mesh nodes

global bacterialEdges;
global bactNodes;
global bactNodeCoeffs;
[p,e,t]= meshToPet(mesh);

% Find the mesh nodes that correspond to each bacterium in the geometry
nBact= length(bacterialEdges);
bactNodes= cell(nBact,1);
% Find bacterial nodes
for b= 1:nBact
  % all nodes for edges on bacterium
  bactNodes{b}= unique( e(1,ismember(e(5,:),bacterialEdges{b})) );
  if isempty(bactNodes{b}), error(['No mesh nodes inside bacterium #',num2str(b)]); end
  fprintf('Nodes in bact#%d: %d\n', b,length(bactNodes{b}));
end

% Calculate dAHL mesh node coefficients
% Each node is part of many edges. Each edge has a length relative to the bacterium's perim
% Each node's relative length is the sum of the relative lengths of each edge it is part of
% divided by the number of nodes for each edge (2)
bactNodeCoeffs= cell(nBact,1);
assert(strcmp(mesh.GeometricOrder,'linear'));
for b= 1:nBact
  nodes= bactNodes{b};
  selEdge= ismember(e(5,:),bacterialEdges{b});
  nodesForeachElt= t(1:3,selEdge);
  nodePerElt= sum( ismember(nodesForeachElt,nodes) );   % number of (selected) nodes foreach elt
  eltRelArea= pdetrg(p,nodesForeachElt)./prod(bactSize);
  nodeRelArea= zeros(length(nodes),1);
  for n= 1:length(nodes)
    nodePartofElt= any(nodesForeachElt == nodes(n));
    nodeRelArea(n)= sum(eltRelArea(nodePartofElt)./nodePerElt(nodePartofElt));
  end
  bactNodeCoeffs{b}= bactRingDensity(b).*nodeRelArea;
end

%ahlProd= dy(6).*bactRingDensity(b).*([1;2;1;2;3;3]./3*0.25);  %<dydt>*<elt/node>/<nodes/elt>*<eltRelArea>
