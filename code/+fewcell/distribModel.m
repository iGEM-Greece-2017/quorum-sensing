% 4th QS model: model a petri dish as a cylinder and solve on 2D due to
% angular symmetry. Algebraically transform the resulting diffusion equation
% to form an equivalent 2d cartesian problem with radial coefficients
% - Help topic: "Heat Distribution in a Circular Cylindrical Rod"
clear params;
%% Parameters
% Units:
% Molarity: nM, Time: min, Length: mm, Quantity: fmol (nM*mm^3)
runID= randi(2147483644);
% time
params.t.tstop= 38*1;   % min
params.t.tstart= 0;
params.t.timePoints= 100;
% coefficients
params.c.c_agar= 1e-3*60;        % [mm^2/min]
params.c.c_cytoplasm= 1e-1*60;   % [mm^2/min]
params.c.d_AHL= 0.01;            % [1/min]
params.c.bactMperm= 100;         % [1/min]
% geometry
params.g.bactCenters= 1e-3*[0,-3];
params.g.bactCenters(2,:)= params.g.bactCenters(1,:)-[4,10.25]*1e-3;
%params.g.bactCenters(3,:)= params.g.bactCenters(1,:)-[-4,10.25]*1e-3;
params.g.bactSize= 1e-3*[1,2];
params.g.domainLim= [10,2];
% mesh
params.m.Hgrad= 1.35;
params.m.HmaxCoeff= 1/15;
% init
params.solve.y0= [1.6;0;0;0;0; 0;0;0];  % [nM]
% solve
params.solve.AbsTol_y= [1e-5,1e-6,1e-5,1e-3,1e-3,  1e-6,1e-6,1e-6].*1e-1;
params.solve.AbsTol= 1e-2;    % for diffusion nodes
params.solve.RelTol= 1e-4;
params.solve.FeatureSize= min(params.g.bactSize)/2;
% viz
params.viz.showMesh= true;
params.viz.zoominFactor= [40,8];
params.viz.interpResolution= 120;
params.viz.integrateAbstol= 1;
params.viz.timePoints= floor(params.t.timePoints/2);
params.viz.dynamicScaling= false;
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
