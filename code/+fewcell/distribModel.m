% 4th QS model: model a petri dish as a cylinder and solve on 2D due to
% angular symmetry. Algebraically transform the resulting diffusion equation
% to form an equivalent 2d cartesian problem with radial coefficients
% - Help topic: "Heat Distribution in a Circular Cylindrical Rod"
%% Parameters
% time
params.t.tstop= 60*1;  % sec
params.t.tstart= 0;
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
params.m.Hgrad= 1.25;
params.m.HmaxCoeff= 1/12;
% solve
params.solve.AbsTol= 1e-10; % for diffusion nodes
params.solve.RelTol= 1e-7;
params.solve.FeatureSize= min(params.g.bactSize);
% viz
params.viz.zoominFactor= [4,2];
params.viz.interpResolution= 80;
params.viz.timePoints= 120;
params.viz.dynamicScaling= true;

%% Solve
[model,tlist,domainVolume]= fewcell.problemSetup(params,true);
fprintf('--> Solving...\n');
result= fewcell.solveProblem(model,tlist,params.solve);

%% Plot solution
% Prepare solution interpolation
fprintf('--> Interpolating solution...\n');
[AHLDistrib,x,y,interp_t,totalAHL]= fewcell.util.interpolateIntegrateAHL(model,result,params);
% Plot
fprintf('--> [paused] Press key to plot the solution\n');
pause;
fprintf('--> Plotting solution...\n');
global distribAHL_interp_graphics;
distribAHL_interp_graphics= []; distribAHL_interp_graphics.first= true;
for t= 1:length(interp_t)
  fewcell.viz.plot(AHLDistrib,x,y,interp_t,totalAHL,3,t);
end
% Printed results
fprintf('Final total [AHL] at t=%.1fsec: %.3g\n\n', interp_t(end), totalAHL(end));
