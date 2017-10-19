% 4th QS model: model a petri dish as a cylinder and solve on 2D due to
% angular symmetry. Algebraically transform the resulting diffusion equation
% to form an equivalent 2d cartesian problem with radial coefficients
% - Help topic: "Heat Distribution in a Circular Cylindrical Rod"
matlab_start;
clear params; clear global;
%% Parameters
% Units:
% Molarity: nM, Time: min, Length: mm, Quantity: fmol (nM*mm^3)
params.runID= randi(2147483644);
global enableSinglecellEq; enableSinglecellEq= true;  % false: debugging only
global enableGraphics; enableGraphics= true;
% time
params.t.tstop= 60*16;   % min
params.t.tstart= 0;
params.t.timePoints= 600;
% coefficients
%params.c.c_agar= 7.1e-5*60;                 % [mm^2/min]
params.c.c_agar= 25e-5*60;                 % [mm^2/min]
params.c.c_cytoplasm= 30e-5*60;            % [mm^2/min]
params.c.d_AHL= 7e-5;                      % [1/min]
%params.c.d_AHL= 0.001;

% geometry
params.g.bactSize= 1e-3*[1,2.164];
params.g.bactCenter0= 1e-3*[300,-1.082];
params.g.nRings= 20; params.g.nLayers= 2;
params.g.ringDist= 5;   % must be an odd number
params.g.layerSeparation= 1;
%params.g.domainLim= [17,5.51];       % small disk
%params.g.domainLim= [1.7, .551];     % xs
params.g.domainLim= [1.4, .45];     % xs
%params.g.domainLim= [.4,.1];         % tiny
%params.g.domainLim= [20,20]*1e-3;    % tiny

% Membrane r: 12.5mm -> (30oC)
% - Full nutrient agar: 10^10.1 cfu/membrane, growth.r=1.71
% *- 1/5 nutrient: 10^9.45 cfu/membrane, growth.r=1.5
% - 1/25 nutrient: 10^8.9 cfu/membrane, growth.r=1.5
% - no nutrients: 10^7.75 cfu/membrane, growth.r=.98

% growth
params.growth.on= false;    % enable/disable growth

% mesh
params.m.Hgrad= 1.5;
params.m.HmaxCoeff= 1/12;
% init
params.solve.y0= [1.5347;0;0;0;0; 0;0;0];  % [nM]
% solve
params.solve.AbsTol_y= [1e-3,1e-2,1e-2,1e-3,1e-3,  1e-6,1e-6,1e-4]*1e0;
params.solve.AbsTol= 1e-7;    % for diffusion nodes
params.solve.RelTol= 1e-4;
params.solve.FeatureSize= min(params.g.bactSize)/10;

% viz
params.viz.showMesh= true;
params.viz.showGrowthCurve= true;
params.viz.domLim= [[0;0],params.g.domainLim'./1];
params.viz.interpResolution= 120;
params.viz.integrateAbstol= 1;
params.viz.timePoints= floor(params.t.timePoints/12);
params.viz.dynamicScaling= true;
params.viz.logscaleSinglecell= false;
params.viz.figID= [5,6,7];

%% Solve
[params,tlist,bactRingDensity]= fewcell.initSetup(params);
[params,model,domainVolume]= fewcell.problemSetup(params);
fprintf('--> Solving...\n');
result= fewcell.solveProblem(model,tlist,params);

%% Plot solution
% Prepare solution interpolatms.g.ringDist= 1*1e-3; params.g.layerSeparation= 1.082*1e-3;
fprintf('--> Interpolating solution...\n');
[AHLDistrib,x,y,interp_t,totalAHL]= fewcell.util.interpolateIntegrateAHL(model,result,params);
% Plot
%if params.viz.figID(1), fprintf('--> [paused] Press key to plot the solution\n'); pause; end
fprintf('--> Plotting solution...\n');
fewcell.output(model,result,tlist,AHLDistrib,x,y,interp_t,totalAHL,params.viz);
