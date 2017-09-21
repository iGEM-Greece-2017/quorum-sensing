function [p,model,tlist,domainVolume,bactRingDensity]= problemSetup(p,plotMesh)   %param
% Setup all elements of the pde problem

%% Time
%ntime= (p.t.tstop - p.t.tstart)+1;
tlist= linspace(p.t.tstart, p.t.tstop, p.t.timePoints);
%% Geometry
% Develop bactCenters from metaparameters
p.g.bactCenters= zeros(p.g.nRings*p.g.nLayers,2);
p.g.bactCenters(1,:)= p.g.bactCenter0;
prevLayerStart= 1;
for layer= 1:p.g.nLayers
  if layer>1
    p.g.bactCenters((layer-1)*p.g.nRings+1,:)= ...
      p.g.bactCenters(prevLayerStart,:)-[0,p.g.layerSeparation+p.g.bactSize(2)];
  end
  for ring= (layer-1)*p.g.nRings+2:layer*p.g.nRings
    p.g.bactCenters(ring,:)= ...
      p.g.bactCenters(ring-1,:)+[p.g.ringDist+p.g.bactSize(1),0];
  end
  prevLayerStart= (layer-1)*p.g.nRings+1;
end

nBact= size(p.g.bactCenters,1);
geometryFun= @(varargin)fewcell.bactAgarGeom(...
                p.g.bactCenters,p.g.bactSize, p.g.domainLim, varargin);
domainVolume= pi*p.g.domainLim(1).^2 * p.g.domainLim(2);  % domain is actually a cylinder
bactVolume= pi/4*p.g.bactSize(1).^2 * p.g.bactSize(2);

bactRingDensity= fewcell.util.bactRingDensity(p.g.bactCenters(:,1),p.g.bactSize, p.g.lateralSpacing);
totalBacteria= round(sum(bactRingDensity));
fprintf('Total bacteria: %d\n', totalBacteria);
%% PDE
model= createpde(1);
geometryFromEdges(model,geometryFun);
%% Boundaries
applyBoundaryCondition(model, 'neumann', 'Edge',1,'g',0,'q',0);
%applyBoundaryCondition(model, 'dirichlet', 'Edge',1,'u',0);
applyBoundaryCondition(model, 'neumann', 'Edge',2:4,'g',0,'q',0);
%% Coefficients
% The del operator is expressed in cylindrical coordinates. du/dθ=0 due to symmetry, so:
%    u' - del*( c del(u)) +  au =  f
% becomes, after substituting r<-x,
%   xu' - del*(xc del(u)) + xau = xf
% c: diffusion, a: destruction rate
dCoeff= @(r,s)r.x;
cCoeff= @(r,s)p.c.c_agar*r.x;
aCoeff= @(r,s)p.c.d_AHL*r.x;
specifyCoefficients(model,'Face',1,'m',0,'d',dCoeff,'c',cCoeff,'a',aCoeff,'f',0);
for b=1:nBact
  %bactProd= @(r,s)(10*s.u+(1e1*(s.time<2)));
  %bactProd= @(r,s)(1*(s.time==tstart));
  %bactProd= @(r,s)1e4*(r.x+s.u);
  % Spatial coefficient for production through singlecell eqs
  %bactProd= @(r,s)2*pi*r.x*p.g.bactSize(1)*p.g.bactSize(2)./bactVolume;
  bactProd= @(r,s)r.x;
  dCoeff= @(r,s)r.x;
  cCoeff= @(r,s)p.c.c_cytoplasm*r.x;
  specifyCoefficients(model,'Face',b+1,...
    'm',0,'d',dCoeff,'c',cCoeff,'a',0,'f',bactProd);
end
%% Initial conditions
setInitialConditions(model,0);
%% Mesh
generateMesh(model,'MesherVersion','R2013a', 'Jiggle','minimum','JiggleIter',50,...
  'Hgrad',p.m.Hgrad, 'Hmax',p.g.domainLim(1) * p.m.HmaxCoeff);
totalMeshNodes= size(model.Mesh.Nodes,2);
fprintf('Total mesh nodes: %d\n', totalMeshNodes);

% Plot mesh
if plotMesh
  figure(2); clf;
  pdeplot(model,'NodeLabels','off');
  %pdegplot(geometryFun, 'EdgeLabels','on','FaceLabels','on');
  axis tight;
  meshPlot= gca;
  meshPlot.XLim= [0 p.viz.domLim(1)];
  meshPlot.YLim= [-p.viz.domLim(2) 0];
  drawnow;
end
