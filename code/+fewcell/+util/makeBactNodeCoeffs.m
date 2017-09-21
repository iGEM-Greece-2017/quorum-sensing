function makeBactNodeCoeffs(mesh,nBact)
% Find the mesh nodes that correspond to each bacterium in the geometry and calculate the 
% coefficient that translates the singlecell dAHL into a change for the corresponding mesh nodes

global bactNodes;
[p,~,t]= meshToPet(mesh);

% Find the mesh nodes that correspond to each bacterium in the geometry
if nargin<4, nBact= length(unique(t(4,:)))-1; end
bactNodes= cell(nBact,1);
% Find bacterial nodes
for b= 1:nBact
  bactNodes{b}= unique( t(1:3,t(4,:)==b+1) );         % all nodes for elements in bacterium
  if isempty(bactNodes{b}), error(['No mesh nodes on bacterium #',num2str(b)]); end
  fprintf('Nodes in bact#%d: %d\n', b,length(bactNodes{b}));
end
