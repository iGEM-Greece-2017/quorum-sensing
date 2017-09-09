function [p,model,tlist,domainVolume]= problemSetup(p,plotMesh)   %param
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
%% PDE
model= createpde(1);
geometryFromEdges(model,geometryFun);
%% Boundaries
global bacterialEdges;
applyBoundaryCondition(model, 'neumann', 'Edge',1,'g',0,'q',0);
applyBoundaryCondition(model, 'neumann', 'Edge',2:4,'g',0,'q',0);
for b=1:nBact
  bacterialEdges{b}= (1:4)+4+(b-1)*4;
  %membranePerm= @(r,s)p.c.bactMperm*r.x;
  %applyBoundaryCondition(model,'neumann', 'Edge',(1:4)+4+(b-1)*4,'g',0,'q',membranePerm);
  membU= @(r,s)(1+s.u.^0.1).*r.x;
  applyBoundaryCondition(model, 'dirichlet', 'Edge',bacterialEdges{b},'h',1,'r',membU);
end
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
  pdeplot(model,'NodeLabels','on');
  %pdegplot(geometryFun, 'EdgeLabels','on','FaceLabels','on');
  axis tight;
  meshPlot= gca;
  meshPlot.XLim= [0 0.02];
  meshPlot.YLim= [-0.01 0];
  drawnow;
end
