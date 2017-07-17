%% Parameters
y0= [1.6E-7;0;0;0;0;0;0;0;0]*1e-9;   % singlecell initial conditions
tstop= 60*5000;  % sec

bactSize= [1.5;1]*1e-6;
bactCenters= [0,0]*1e-6;
domainLim= 1e-4;
zoominFactor= 1;    % When visualizing, focus on a small part of the domain

interpResolution= 80;
timeSubsampling= tstop/80;
dynamicScaling= true;

%% Time discretization
ntime= tstop+1;
tlist= linspace(30,tstop,ntime);
%% Geometry
geometryFun= @(varargin)fewcell.topview.bactDomainGeom(bactCenters,bactSize,domainLim,varargin);
domainArea= (2*domainLim + domainLim) * sin(pi/3)*domainLim;  % hexagon:= 2*trapezoids
%% PDE
numPde= 1;
model= createpde(numPde);
geometryFromEdges(model,geometryFun);
%{
% Plot pde domain
figure(1);
pdegplot(model);%,'EdgeLabels','on','FaceLabels','on');
drawnow;
axis equal;
%}

%% Boundaries
nBact= size(bactCenters,1);
applyBoundaryCondition(model,'neumann','Edge',1:6,'g',0,'q',0);
for b= 1:nBact
  applyBoundaryCondition(model,'neumann','Edge',(1:2)+6+(b-1)*2,'g',0,'q',100);
end
%% PDE coeffs
%prodFun= @(~,st)(exp(-st.u*1e10)-0.85)*20e-6;
%bacteriaTable= fewcell.topview.BacteriaTable(1,y0,tlist(1));
% d: d*u', c: -div(c*grad(u)), f: source
specifyCoefficients(model, 'Face',1, 'm',0, 'd',1, 'c',1e-9, 'a',1/3, 'f',0);
for b= 1:nBact
  % Dummy AHL production. Helps identify the nodes of the bacterium
  specifyCoefficients(model, 'Face',b+1, 'm',0, 'd',1, 'c',1e-7, 'a',0, 'f',@(r,s)b*(s.u+1));
end
%% Initial conditions
setInitialConditions(model,0);
%% Mesh
mesh= generateMesh(model,'MesherVersion','R2013a', 'Jiggle','minimum', 'JiggleIter',25,...
  'GeometricOrder','linear', 'Hgrad',1.9,'Hmax',domainLim /6);
totalMeshNodes= size(model.Mesh.Nodes,2);
fprintf('Total mesh nodes: %d\n', totalMeshNodes);
%{
% Plot mesh
figure(2); clf;
pdeplot(model);%,'NodeLabels','on');
drawnow;
axis equal;
%}

%% Solve
fprintf('--> Solving...\n');
tic;
model.SolverOptions.AbsoluteTolerance= 1e-25;
model.SolverOptions.RelativeTolerance= 1e-12;
model.SolverOptions.ResidualTolerance= 1e-12;
model.SolverOptions.MaxIterations= 50;
model.SolverOptions.MinStep= min(bactSize)/10;
model.SolverOptions.ReportStatistics= 'on';
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
  tic;
  fewcell.viz.distribAHL_interp(AHLDistrib,x,y,totalAHL,interpTimes,...
    domainLim/zoominFactor,dynamicScaling,3,t);
  % Control framerate
  pauseTime= 1/62 - toc;
  if pauseTime>1e-4, pause(pauseTime); end
end

% Printed results
% ref: http://book.bionumbers.org/what-is-the-concentration-of-bacterial-cells-in-a-saturated-culture/
bacterialDensity= 1/(domainArea*1e4); % [bacteria/cm^2]
fprintf('Bacterial density: %.2g bacteria/cm^2\t %.2g bacteria/mL (stacked tightly)\n', bacterialDensity, bacterialDensity*1e-2/(2*min(bactSize)));
fprintf('Final total [AHL] at t=%.1fsec: %.3g\n\n', interpTimes(end), totalAHL(end));
