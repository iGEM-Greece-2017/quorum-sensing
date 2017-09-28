%% Estimate an equivalent diffusion constant to compensate for thet membrane permeability
matlab_start;
global enableSinglecellEq; enableSinglecellEq=true;

% ** DEBUG
global prodDist;
prodDist.t= []; prodDist.F= []; prodDist.femF= [];
global dyLog;
dyLog.y= []; dyLog.t= [];

%% Parameters
tlist= linspace(0, 60*16, 400);
params.g.bactSize= [1,2.164]*1e-3;
params.g.domainLim= [1,2]*1e-2;
params.c.high= 1e-6;
params.c.dAHL= 0;

params.viz.domLim= params.g.bactSize/2; params.viz.zoominFactor= params.g.domainLim./params.viz.domLim;
params.viz.interpResolution= 150; params.viz.timePoints= floor(length(tlist)/4);
params.viz.figID= [0,3,4]; params.viz.integrateAbstol= 1;
params.viz.dynamicScaling= true; params.viz.logscaleSinglecell= false;

%% Geometry: 3 rectangles
x= [1,0;0,1;0,1;1,0]; y= [1,0;1,0;0,1;0,1];
shapes= [3;4; x*[0;params.g.bactSize(1)/2]; y*[0;-params.g.bactSize(2)/2]]; %bactOnly
names= char('bact'); names= names'; setf= 'bact';

[g,~]= decsg(shapes,setf,names);
domainVol= (params.g.domainLim(1)/2)^2*pi*params.g.domainLim(2);
fprintf('Domain volume: %e\n', domainVol);
%% BC
model= createpde(1);
geometryFromEdges(model,g);
%figure(6); pdegplot(model,'EdgeLabels','on','FaceLabels','on');

applyBoundaryCondition(model, 'neumann', 'Edge',1:4,'g',0,'q',0);

%% Coeffs
% Test relationship between nodal F and analytical production (linear)
prodCoeff= @(r,s)r.x;    % calibrates nodal coefficients
specifyCoefficients(model,'m',0,'d',@(r,s)r.x,'c',@(r,s)r.x*params.c.high,'a',@(r,s)r.x*params.c.dAHL,'f', prodCoeff);    %bact
%specifyCoefficients(model,'m',0,'d',1,'c',params.c.high,'a',params.c.dAHL,'f', prodCoeff);    %bact

setInitialConditions(model,0);

generateMesh(model,'MesherVersion','R2013a','GeometricOrder','linear', 'Jiggle','minimum','JiggleIter',50,...
  'Hgrad',1.3, 'Hmax',params.g.domainLim(1)/40);
totalMeshNodes= size(model.Mesh.Nodes,2);
fprintf('Total mesh nodes: %d\n', totalMeshNodes);
%figure(6); clf; pdeplot(model, 'NodeLabels','on'); drawnow;

%% Solve
global bactNodes; bactNodes= true(totalMeshNodes,1);

bactInterior= params.g.bactSize;%-2*params.g.membW;
bactArea= prod(bactInterior)/4;
bactVol= (params.g.bactSize(1)/2)^2*pi*params.g.bactSize(2);
global bactWallDilution; bactWallDilution= 1;%/(prod(params.g.bactSize)/prod(bactInterior));

model.SolverOptions.AbsoluteTolerance= 1e-5;
model.SolverOptions.RelativeTolerance= 1e-4;
model.SolverOptions.ResidualTolerance= 1e-4;
model.SolverOptions.MaxIterations= 40;
model.SolverOptions.MinStep= params.g.bactSize(1)/20;
model.SolverOptions.ReportStatistics= 'on';

global solveInternalParams;
solveInternalParams.y0= [1.5347;0;0;0;0; 0;0;0];
solveInternalParams.AbsTol_y= [1e-5,1e-6,1e-5,1e-3,1e-3,  1e-6,1e-6,1e-7]*1e-1;
tic;
result= solvepde(model,tlist);
toc;
cAHL= result.NodalSolution(5,:)';
global yyResults;
for b= 1:length(yyResults)
  yyResults{b}(:,[6,10])= [cAHL,zeros(length(tlist),1)];
  yyResults{b}= [yyResults{b},ones(length(tlist),1)];
end

%% Display
fprintf('--> Interpolating solution...\n');
[AHLDistrib,x,y,interp_t,totalAHL]= ...
  fewcell.util.interpolateIntegrateAHL(model,result,params);
% Plot
fprintf('--> Plotting solution...\n');
fewcell.output(model,result,tlist,...
               AHLDistrib,x,y,interp_t,totalAHL,params.viz);

%figure(5); plot(tlist/60,cAHL);

