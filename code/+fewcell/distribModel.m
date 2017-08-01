% 4th QS model: model a petri dish as a cylinder and solve on 2D due to
% angular symmetry. Algebraically transform the resulting diffusion equation
% to form an equivalent 2d cartesian problem with radial coefficients
% - Help topic: "Heat Distribution in a Circular Cylindrical Rod"
%% Parameters
% time
params.t.tstop= 60*60*1;  % sec
params.t.tstart= 0;
params.t.timePoints= 30;
% coefficients
params.c.c_agar= 1e-9;
params.c.c_cytoplasm= 1e-7;
params.c.d_AHL= 0.01/60;
params.c.bactMperm= 100;
% geometry
params.g.bactCenters= [0,-10]*1e-6;
params.g.bactSize= [2,1]*1e-6;
params.g.domainLim= [0.5e-2,2e-3];
% mesh
params.m.Hgrad= 1.3;
params.m.HmaxCoeff= 1/12;
% solve
params.solve.AbsTol= 1e-8; % for diffusion nodes
params.solve.RelTol= 1e-8;
params.solve.FeatureSize= min(params.g.bactSize)/2;
% viz
params.viz.zoominFactor= [5,3];
params.viz.interpResolution= 100;
params.viz.timePoints= params.t.timePoints;
params.viz.dynamicScaling= true;

%% Solve
[model,tlist,domainVolume]= fewcell.problemSetup(params,true);
params.viz.timePoints= length(tlist);
fprintf('--> Solving...\n');
result= fewcell.solveProblem(model,tlist,params.solve);

%% Plot solution
% Prepare solution interpolation
fprintf('--> Interpolating solution...\n');
[AHLDistrib,x,y,interp_t,totalAHL]= fewcell.util.interpolateIntegrateAHL(model,result,params);
% Plot
fprintf('--> [paused] Press key to plot the solution\n'); pause;
fprintf('--> Plotting solution...\n');
fewcell.output(AHLDistrib,x,y,interp_t,totalAHL,false);
