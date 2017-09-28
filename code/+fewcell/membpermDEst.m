%% Estimate an equivalent diffusion constant to compensate for thet membrane permeability
matlab_start;
global enableSinglecellEq; enableSinglecellEq=true;

% ** DEBUG
global prodDist;
prodDist.t= []; prodDist.F= []; prodDist.femF= [];

%% Parameters
tlist= linspace(0, 60*16, 400);
params.g.membW= 0.2e-3;
params.g.bactSize= [1,2.164]*1e-3;
%params.g.bactSize= [.8,3.38125]*1e-3;    % constant volume (only change aspect ratio
params.g.domainLim= params.g.bactSize*5;
params.c.high= 1e-3;
params.c.wall= 1e-6;
params.c.dAHL= .001;

params.viz.domLim= params.g.domainLim/2; params.viz.zoominFactor= params.g.domainLim./params.viz.domLim;
params.viz.interpResolution= 150; params.viz.timePoints= floor(length(tlist)/6);
params.viz.figID= [2,3,4]; params.viz.integrateAbstol= 1e-1;
params.viz.dynamicScaling= true; params.viz.logscaleSinglecell= false;

%% Geometry: 3 rectangles
x= [1,0;0,1;0,1;1,0]; y= [1,0;1,0;0,1;0,1];
shapes(:,1)= [3;4; x*[0;params.g.domainLim(1)/2]; y*[0;-params.g.domainLim(2)/2]];  %dom
shapes(:,2)= [3;4; x*[0;params.g.bactSize(1)/2]; y*[0;-params.g.bactSize(2)/2]];  %bact
shapes(:,3)= [3;4; x*[0;params.g.bactSize(1)/2+params.g.membW]; y*[0;-params.g.membW-params.g.bactSize(2)/2]]; %memb
names= char('dom','bact','memb'); names= names';
setf= 'bact+dom+memb';

[g,~]= decsg(shapes,setf,names);
bactVol= (params.g.bactSize(1)/2)^2*pi*params.g.bactSize(2);
domainVol= (params.g.domainLim(1)/2)^2*pi*params.g.domainLim(2);
fprintf('Bact volume: %.4e\tDomain volume: %.4e\n', bactVol, domainVol);
%% BC
model= createpde(1);
geometryFromEdges(model,g);
%figure(6); pdegplot(model,'EdgeLabels','on','FaceLabels','on');

applyBoundaryCondition(model, 'neumann', 'Edge',1:12,'g',0,'q',0);

%% Coeffs
dCoeff= @(r,s)r.x;
aCoeff= @(r,s)params.c.dAHL*r.x;
% Test relationship between nodal F and analytical production (linear)
%prodCoeff= @(r,s) (1e2*(s.time<1)+s.time*1e2*(s.time>=1))*ones(size(s.u));
prodCoeff= @(r,s)r.x;    % calibrates nodal coefficients
specifyCoefficients(model,'Face',1,'m',0,'d',dCoeff,'c',@(r,s)r.x*params.c.high,'a',aCoeff,'f',0);     %dom
specifyCoefficients(model,'Face',2,'m',0,'d',dCoeff,'c',@(r,s)r.x*params.c.high,'a',0,'f', prodCoeff);    %bact
specifyCoefficients(model,'Face',3,'m',0,'d',dCoeff,'c',@(r,s)r.x*params.c.wall,'a',aCoeff,'f',0);     %memb

setInitialConditions(model,0);

generateMesh(model,'MesherVersion','R2013a','GeometricOrder','linear', 'Jiggle','minimum','JiggleIter',50,...
  'Hgrad',1.01, 'Hmax',params.g.domainLim(1)/25);
totalMeshNodes= size(model.Mesh.Nodes,2);
fprintf('Total mesh nodes: %d\n', totalMeshNodes);
figure(6); clf; pdeplot(model, 'NodeLabels','off'); drawnow;

%% Solve
fewcell.util.makeBactNodeCoeffs(model.Mesh,1);

bactArea= prod(params.g.bactSize)/4;
bactVol= (params.g.bactSize(1)/2)^2*pi*params.g.bactSize(2);

model.SolverOptions.AbsoluteTolerance= 1e-8;
model.SolverOptions.RelativeTolerance= 1e-5;
model.SolverOptions.ResidualTolerance= 1e-5;
model.SolverOptions.MaxIterations= 40;
model.SolverOptions.MinStep= params.g.membW/20;
model.SolverOptions.ReportStatistics= 'on';

global solveInternalParams;
solveInternalParams.y0= [1.5347;0;0;0;0; 0;0;0];
solveInternalParams.AbsTol_y= [1e-5,1e-6,1e-5,1e-3,1e-3,  1e-6,1e-6,1e-7]*1e-1;
tic;
result= solvepde(model,tlist);
toc;
cAHL= result.NodalSolution([10,1],:)';
global yyResults;
for b= 1:length(yyResults)
  yyResults{b}(:,[6,10])= cAHL;
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

