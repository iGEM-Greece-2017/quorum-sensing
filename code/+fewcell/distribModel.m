% 4th QS model: model a petri dish as a cylinder and solve on 2D due to
% angular symmetry. Algebraically transform the resulting diffusion equation
% to form an equivalent 2d cartesian problem with radial coefficients
% - Help topic: "Heat Distribution in a Circular Cylindrical Rod"
%% Parameters
runID= randi(2147483644);
% time
params.t.tstop= 50*1*1;  % sec
params.t.tstart= 0;
params.t.timePoints= 10;
% coefficients
params.c.c_agar= 1e-9;
params.c.c_cytoplasm= 1e-7;
params.c.d_AHL= 0.01/60;
params.c.bactMperm= 100;
% geometry
params.g.bactCenters= 1e-6*[0,-2; 0,-5];
params.g.bactSize= 1e-6*[1,2];
params.g.domainLim= [0.5e-2,2e-3] ./100;
% mesh
params.m.Hgrad= 1.05;
params.m.HmaxCoeff= 1/20;
% init
params.solve.y0= [1.6E-9;0;0;0;0; 0;0;0];
% solve
%params.solve.AbsTol_y= [1e-10,5e-12,1e-11,1e-9,1e-8, 1e-13,1e-16,1e-18]*1e0;
params.solve.AbsTol_y= [1e-10,5e-12,1e-11,1e-10,1e-11, 1e-13,1e-26,1e-28]*1e-2;
params.solve.AbsTol= 1e-17;  % for diffusion nodes
params.solve.RelTol= 1e-9;
params.solve.FeatureSize= min(params.g.bactSize)/10;
% viz
params.viz.showMesh= true;
params.viz.zoominFactor= [1,1];
params.viz.interpResolution= 130;
params.viz.timePoints= floor(params.t.timePoints/2);
params.viz.dynamicScaling= true;
%params.viz.figID= runID+[1,2];
params.viz.figID= [3,4];

%% Solve
[model,tlist,domainVolume]= fewcell.problemSetup(params,params.viz.showMesh);
fprintf('--> Solving...\n');
result= fewcell.solveProblem(model,tlist,params.solve);

%% Plot solution
% Prepare solution interpolation
fprintf('--> Interpolating solution...\n');
[AHLDistrib,x,y,interp_t,totalAHL]= fewcell.util.interpolateIntegrateAHL(model,result,params);
% Plot
fprintf('--> [paused] Press key to plot the solution\n'); pause;
fprintf('--> Plotting solution...\n');
fewcell.output(AHLDistrib,x,y,interp_t,totalAHL,params.viz.figID,params.viz.dynamicScaling);

global yyResults;
singlecell.viz.plot(tlist, yyResults{1}, 1,params.viz.figID(2));
