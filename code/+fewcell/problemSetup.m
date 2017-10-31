function [p,model,domainVolume]= problemSetup(p,initConditionsFun)   %param
% Setup all elements of the pde problem
global enableGraphics
%% Geometry
domainVolume= pi*p.g.domainLim(1).^2 * p.g.domainLim(2);  % domain is actually a cylinder
bactVolume= pi/4*p.g.bactSize(1).^2 * p.g.bactSize(2);
% First bactCenter
p.g.bactCenters= zeros(p.g.nRings*p.g.nLayers,2);
bactAligners= zeros([(p.g.ringDist-1)/2,p.g.nRings*p.g.nLayers,2]);
p.g.bactCenters(1,:)= p.g.bactCenter0;
for s= 1:(p.g.ringDist-1)/2
  bactAligners(s,1,:)= ...
    p.g.bactCenters(1,:)+[s*2*p.g.bactSize(1),0];
end
% Develop the rest of the bactCenters from metaparameters
ringDist= p.g.ringDist*p.g.bactSize(1);
assert(mod(p.g.ringDist,2)==1, 'Ring distance must be an odd number');
layerSeparation= p.g.layerSeparation*p.g.bactSize(2);
prevLayerStart= 1;
for layer= 1:p.g.nLayers
  if layer>1
    p.g.bactCenters((layer-1)*p.g.nRings+1,:)= ...
      p.g.bactCenters(prevLayerStart,:)-[0,layerSeparation+p.g.bactSize(2)];
    for s= 1:(p.g.ringDist-1)/2
      bactAligners(s,(layer-1)*p.g.nRings+1,:)= ...
        p.g.bactCenters((layer-1)*p.g.nRings+1,:)+[s*2*p.g.bactSize(1),0];
    end
  end
  for ring= (layer-1)*p.g.nRings+2:layer*p.g.nRings
    p.g.bactCenters(ring,:)= ...
      p.g.bactCenters(ring-1,:)+[ringDist+p.g.bactSize(1),0];
    for s= 1:(p.g.ringDist-1)/2
      bactAligners(s,ring,:)= ...
        p.g.bactCenters(ring,:)+[s*2*p.g.bactSize(1),0];
    end
  end
  prevLayerStart= (layer-1)*p.g.nRings+1;
end
bactAligners= reshape(bactAligners,[],2);
nBact= size(p.g.bactCenters,1);

% Create geometry description
x= [1,0;0,1;0,1;1,0]; y= [1,0;1,0;0,1;0,1];
names= cell(2,nBact+1);
names(:,1)= {'domain', '+'};
shapes= zeros(10,nBact+1);
shapes(:,1)= [3;4; x*[0;p.g.domainLim(1)]; y*[0;-p.g.domainLim(2)]];  %domain
% Bacterial geometry
for b= 1:nBact
  bsize= p.g.bactSize;
  bULcorner= p.g.bactCenters(b,:) - [bsize(1)/2,-bsize(2)/2];   % upper left corner
  shapes(:,b+1)= [3;4; x*[bULcorner(1); bULcorner(1)+bsize(1)]; y*[bULcorner(2); bULcorner(2)-bsize(2)]];
  names(:,b+1)= {['bact',num2str(b)]; '+'};
end
for s= 1:size(bactAligners,1)
  bsize= p.g.bactSize;
  bULcorner= bactAligners(s,:) - [bsize(1)/2,-bsize(2)/2];   % upper left corner
  shapes(:,nBact+s+1)= [3;4; x*[bULcorner(1); bULcorner(1)+bsize(1)]; y*[bULcorner(2); bULcorner(2)-bsize(2)]];
  names(:,nBact+s+1)= {['algn',num2str(s)]; '+'};
end

% Create the description
names{2,end}= '';
setf= [names{:}];
names= char(names{1,:}); names= names';
[geometryDescription,~]= decsg(shapes,setf,names);

%% PDE
model= createpde(1);
geometryFromEdges(model,geometryDescription);
%% Mesh
generateMesh(model,'MesherVersion','R2013a','GeometricOrder','linear', 'Jiggle','minimum','JiggleIter',50,...
  'Hgrad',p.m.Hgrad, 'Hmax',p.g.domainLim(1) * p.m.HmaxCoeff);
totalMeshNodes= size(model.Mesh.Nodes,2);
fprintf('Mesh nodes: %d\tEquations: %d\n', totalMeshNodes, totalMeshNodes + nBact*8);

% Determine bacterial nodes and faces
bactSubdomain= fewcell.util.makeBactNodeCoeffs(model.Mesh,nBact,p.g);

% Plot mesh
if p.viz.showMesh && enableGraphics
  try set(groot,'CurrentFigure',2); catch, figure(2); end
  clf;
  pdeplot(model,'NodeLabels','off');
  %pdegplot(model, 'EdgeLabels','on','FaceLabels','on');
  axis tight;
  meshPlot= gca;
  meshPlot.XLim= p.viz.domLim(1,:);
  meshPlot.YLim= -flip(p.viz.domLim(2,:));
  drawnow;
end
%% Boundaries
numEdges= model.Geometry.NumEdges;
applyBoundaryCondition(model, 'neumann', 'Edge',1:numEdges,'g',0,'q',0);
%% Coefficients
% The del operator is expressed in cylindrical coordinates. du/dθ=0 due to symmetry, so:
%    u' - del*( c del(u)) +  au =  f
% becomes, after substituting r<-x,
%   xu' - del*(xc del(u)) + xau = xf
% c: diffusion, a: destruction rate
dCoeff= @(r,s)r.x;
cCoeff= @(r,s)p.c.c_agar*r.x;
aCoeff= @(r,s)p.c.d_AHL*r.x;
specifyCoefficients(model,'Face',1,'m',0,'d',dCoeff,'c',cCoeff,'a',aCoeff,'f',0);   % Actual domain face
cCoeff= @(r,s)p.c.c_cytoplasm*r.x;
domainSubfaces= setdiff((2:model.Geometry.NumFaces), bactSubdomain);
specifyCoefficients(model,'Face',domainSubfaces,'m',0,'d',dCoeff,'c',cCoeff,'a',0,'f',0);  % pseudo-bacteria
for b=1:nBact
  %bactProd= @(r,s)(10*s.u+(1e1*(s.time<2)));
  %bactProd= @(r,s)(1*(s.time==tstart));
  %bactProd= @(r,s)1e4*(r.x+s.u);
  % Spatial coefficient for production through singlecell eqs
  %bactProd= @(r,s)2*pi*r.x*p.g.bactSize(1)*p.g.bactSize(2)./bactVolume;
  bactProd= @(r,s)r.x;
  cCoeff= @(r,s)p.c.c_cytoplasm*r.x;
  specifyCoefficients(model,'Face',bactSubdomain(b),...
    'm',0,'d',dCoeff,'c',cCoeff,'a',0,'f',bactProd);
end
%% Initial conditions
if nargin==1, initConditionsFun= 0; end
setInitialConditions(model,initConditionsFun);
