function makeBactNodeCoeffs(mesh,nBact)
% Find the mesh nodes that correspond to each bacterium in the geometry and calculate the 
% coefficient that translates the singlecell dAHL into a change for the corresponding mesh nodes

global bactNodes;
global bactNodesEqulength;
[p,~,t]= meshToPet(mesh);

printPerBact= nBact < 15;
faceIdx= 4;
if strcmp(mesh.GeometricOrder, 'quadratic'), faceIdx= 7; end

%{
% Find the mesh nodes that correspond to each bacterium in the geometry
bactNodes= cell(nBact,1);
% Find bacterial nodes
for b= 1:nBact
  bactNodes{b}= unique( t(1:end-1,t(faceIdx,:)==b+2) );         % all nodes for elements in bacterium
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
  candidates= sparse(unique( t(1:end-1,t(faceIdx,:)==b+1) ), 1, true);
  bactNodes(candidates & ~occupied(1:length(candidates)), b)= true;
  occupied(bactNodes(:,b))= true;
  if ~nnz(bactNodes(:,b)), error(['No mesh nodes on bacterium #',num2str(b)]); end
  if printPerBact, fprintf('Nodes in bact#%d: %d\n', b,nnz(bactNodes(:,b))); end
end
bactNodesEqulength= length(unique(sum(bactNodes,1))) == 1;

if ~printPerBact
  bactNodeN= full(sum(bactNodes,1));
  fprintf('[makeBactNodeCoeffs]: Too many bacteria to list. Printing max/min node counts\n');
  fprintf('Min bact nodes: %d\nMax bact nodes: %d\n', min(bactNodeN), max(bactNodeN));
end
