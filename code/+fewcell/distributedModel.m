%% Parameters
y0= [1.6E-7;0;0;0;0;0;0;0;0];   % singlecell initial conditions
tstop= 60*30;  % sec

bactSize= [2;1]*1e-6;
bactCenters= [0,0];
domainLim= 5e-5;
zoominFactor= 2;    % When visualizing, focus on a small part of the domain

interpResolution= 60;
timeSubsampling= tstop/10;
dynamicScaling= true;

%% Time discretization
ntime= tstop*10+1;
tlist= linspace(0,tstop,ntime);
%% Geometry
geometryFun= @(varargin)fewcell.bactDomainGeom(bactCenters,bactSize,domainLim,varargin);
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
applyBoundaryCondition(model,'neumann','Edge',1:6,'g',0,'q',0);
applyBoundaryCondition(model,'neumann','Edge',7:8,'g',0,'q',100);
%% PDE coeffs
%prodFun= @(~,st)(exp(-st.u*1e10)-0.85)*20e-6;
bacteriaTable= fewcell.BacteriaTable(1,y0,tlist(1));
% d: d*u', c: -div(c*grad(u)), f: source
specifyCoefficients(model, 'Face',1, 'm',0, 'd',1, 'c',1e-9, 'a',0, 'f',0);
specifyCoefficients(model, 'Face',2, 'm',0, 'd',1, 'c',1e-9, 'a',0, 'f',@(r,s)bacteriaTable.AHLProd(r,s));
%% Initial conditions
setInitialConditions(model,0);
%% Mesh
mesh= generateMesh(model,'MesherVersion','preR2013a','GeometricOrder','linear', ...
  'Hgrad',1.4,'Hmax',2*domainLim /12);
totalMeshNodes= size(model.Mesh.Nodes,2)
%
% Plot mesh
figure(2); clf;
pdeplot(model);
drawnow;
axis equal;
%}

%% Solve
fprintf('--> Solving...\n');
tic;
model.SolverOptions.AbsoluteTolerance= 1e-15;
model.SolverOptions.RelativeTolerance= 1e-6;
model.SolverOptions.ResidualTolerance= 1e-7;
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
for t=1:nframes
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
fprintf('Final total [AHL] at t=%.1fsec: %.2g\n\n', interpTimes(end), totalAHL(end));
