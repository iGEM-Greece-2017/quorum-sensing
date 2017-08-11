function [model,tlist,domainVolume]= problemSetup(p,plotMesh)   %param
% Setup all elements of the pde problem

%% Time
%ntime= (p.t.tstop - p.t.tstart)+1;
tlist= linspace(p.t.tstart, p.t.tstop, p.t.timePoints);
%% Geometry
geometryFun= @(varargin)fewcell.bactAgarGeom(...
                p.g.bactCenters,p.g.bactSize, p.g.domainLim, varargin);
domainVolume= pi*p.g.domainLim(1).^2 * p.g.domainLim(2);  % domain is actually a cylinder
%% PDE
model= createpde(1);
geometryFromEdges(model,geometryFun);
%% Boundaries
nBact= size(p.g.bactCenters,1);
applyBoundaryCondition(model, 'neumann', 'Edge',1,'g',0,'q',0);
%applyBoundaryCondition(model, 'dirichlet', 'Edge',2:4,'u',0);
applyBoundaryCondition(model, 'neumann', 'Edge',2:4,'g',0,'q',0);
for b=1:nBact
  applyBoundaryCondition(model,'neumann', 'Edge',(1:4)+4+(b-1)*4,'g',0,'q',p.c.bactMperm);
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
for b=1:nBact
  %bactProd= @(r,s)(10*s.u+(1e1*(s.time<2)));
  %bactProd= @(r,s)1e-12*(s.u+r.x);
  bactProd= 0;
  %bactProd= @(r,s)(1*(s.time==tstart));
  dCoeff= @(r,s)r.x;
  cCoeff= @(r,s)p.c.c_cytoplasm*r.x;
  aCoeff= @(r,s)p.c.d_AHL*r.x;
  specifyCoefficients(model,'Face',b+1,...
    'm',0,'d',dCoeff,'c',cCoeff,'a',aCoeff,'f',bactProd);
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
  pdeplot(model);%,'NodeLabels','on');
  %pdegplot(geometryFun, 'EdgeLabels','on','FaceLabels','on');
  axis tight;
  drawnow;
end
