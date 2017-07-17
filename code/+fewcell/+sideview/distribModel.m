% 4th QS model: model a petri dish as a cylinder and solve on 2D due to
% angular symmetry. Algebraically transform the resulting diffusion equation
% to form an equivalent 2d cartesian problem with radial coefficients
% - Help topic: "Heat Distribution in a Circular Cylindrical Rod"
%% Parameters
% time
tstop= 60*10;  % sec
tstart= 0;
% coefficients
c_agar= 1e-9;
c_cytoplasm= 1e-7;
d_AHL= 0.01/60;
% geometry
bactCenters= [0,-2]*1e-6;
bactSize= [2,1]*1e-6;
domainLim= [1e-2,2e-3]/4;
% viz
zoominFactor= 1;
interpResolution= 200;
timeSubsampling= floor((tstop-tstart)/60);
dynamicScaling= true;

%% Time
ntime= (tstop-tstart)+1;
tlist= linspace(tstart, tstop, ntime);
%% Geometry
geometryFun= @(varargin)fewcell.sideview.bactAgarGeom(...
                bactCenters,bactSize, domainLim, varargin);
domainArea= pi*domainLim(1).^2*domainLim(2);  % domain is actually a cylinder
%% PDE
model= createpde(1);
geometryFromEdges(model,geometryFun);
%% Boundaries
nBact= size(bactCenters,1);
applyBoundaryCondition(model, 'neumann', 'Edge',1:4,'g',0,'q',0);
for b=1:nBact
  applyBoundaryCondition(model,'neumann', 'Edge',(1:2)+4+(b-1)*2,'g',0,'q',100);
end
%% Coefficients
% The del operator is expressed in cylindrical coordinates. du/dθ=0 due to symmetry, so:
%    u' - del*( c del(u)) +  au =  f
% becomes, after substituting r<-x,
%   xu' - del*(xc del(u)) + xau = xf
% c: diffusion, a: destruction rate
dCoeff= @(r,s)r.x;
cCoeff= @(r,s)c_agar*r.x;
aCoeff= @(r,s)d_AHL*r.x;
specifyCoefficients(model,'Face',1,'m',0,'d',dCoeff,'c',cCoeff,'a',aCoeff,'f',0);
for b=1:nBact
  bactProd= @(r,s)(1e-2*s.u+(1e-10*(s.time==tstart)));
  specifyCoefficients(model,'Face',b+1,'m',0,'d',1,'c',c_cytoplasm,'a',d_AHL,'f',bactProd);
end
%% Initial conditions
setInitialConditions(model,0);
%% Mesh
generateMesh(model,'MesherVersion','R2013a','Jiggle','minimum','JiggleIter',30,...
  'Hgrad',1.5,'Hmax',domainLim(1)/12);
totalMeshNodes= size(model.Mesh.Nodes,2);
fprintf('Total mesh nodes: %d\n', totalMeshNodes);
%
% Plot mesh
figure(2); clf;
pdeplot(model);%,'NodeLabels','on');
axis tight;
drawnow;
%}

%% Solve
fprintf('--> Solving...\n');
model.SolverOptions.AbsoluteTolerance= 1e-12;
model.SolverOptions.RelativeTolerance= 1e-5;
model.SolverOptions.ResidualTolerance= 1e-5;
model.SolverOptions.MaxIterations= 50;
model.SolverOptions.MinStep= min(bactSize)/10;
model.SolverOptions.ReportStatistics= 'on';
tic;
result= solvepde(model,tlist);
toc;

%% Plot solution
% Prepare solution interpolation
fprintf('--> Interpolating solution...\n');
tic;
[AHLDistrib,x,y,interpTimes,totalAHL]= ...
  fewcell.util.interpolateIntegrateAHL(model,result,domainLim,zoominFactor, ...
                                       interpResolution, timeSubsampling);
toc;
% Plot
fprintf('--> [paused] Press key to plot the solution\n');
pause;
fprintf('--> Plotting solution...\n');
global distribAHL_interp_graphics;
distribAHL_interp_graphics= []; distribAHL_interp_graphics.first= true;
nframes= length(interpTimes);
for t=10:nframes
  fewcell.viz.distribAHL_interp(AHLDistrib,x,y,totalAHL,interpTimes,...
    domainLim/zoominFactor,dynamicScaling,3,t);
end

% Printed results
fprintf('Final total [AHL] at t=%.1fsec: %.3g\n\n', interpTimes(end), totalAHL(end));
