% 4th QS model: model a petri dish as a cylinder and solve on 2D due to
% angular symmetry. Algebraically transform the resulting diffusion equation
% to form an equivalent 2d cartesian problem with radial coefficients
% - Help topic: "Heat Distribution in a Circular Cylindrical Rod"
clear params;
%% Parameters
% Units:
% Molarity: nM, Time: min, Length: mm, Quantity: fmol (nM*mm^3)
runID= randi(2147483644);
global enableSinglecellEq; enableSinglecellEq= true;
% time
params.t.tstop= 60*50;   % min
params.t.tstart= 0;
params.t.timePoints= 100;
% coefficients
params.c.c_agar= 7.1e-5*60;                % [mm^2/min]
params.c.c_cytoplasm= params.c.c_agar*1e2; % [mm^2/min]
params.c.d_AHL= 0.01;                      % [1/min]
params.c.bactMperm= 100;                   % [1/min]
% geometry
nRings= 1; nLayers= 1; ringDist= 1.5; layerSeparation= 2.5;
bactCenter0= 1e-3*[2,-1.5];
params.g.bactSize= 1e-3*[1,2];
params.g.domainLim= [2.5,0.5];
% mesh
params.m.Hgrad= 1.4;
params.m.HmaxCoeff= 1/20;
% init
params.solve.y0= [1.6;0;0;0;0; 0;0;0];  % [nM]
% solve
params.solve.AbsTol_y= [1e-5,1e-6,1e-5,1e-3,1e-3,  1e-6,1e-6,1e-6].*1e-2;
params.solve.AbsTol= 1e-3;    % for diffusion nodes
params.solve.RelTol= 1e-4;
params.solve.FeatureSize= min(params.g.bactSize)/4;
% viz
params.viz.showMesh= true;
params.viz.zoominFactor= [75,15];
params.viz.interpResolution= 130;
params.viz.integrateAbstol= 1;
params.viz.timePoints= floor(params.t.timePoints/2.6);
params.viz.dynamicScaling= true;
params.viz.figID= [3,4,5];

% develop bactCenters
params.g.bactCenters= zeros(nRings*nLayers,2);
params.g.bactCenters(1,:)= bactCenter0;
prevLayerStart= 1;
for layer= 1:nLayers
  if layer>1
    params.g.bactCenters((layer-1)*nRings+1,:)= params.g.bactCenters(prevLayerStart,:)-[0,layerSeparation]*1e-3;
  end
  for ring= (layer-1)*nRings+2:layer*nRings
    params.g.bactCenters(ring,:)= params.g.bactCenters(ring-1,:)+[ringDist,0]*1e-3;
  end
  prevLayerStart= (layer-1)*nRings+1;
end

%% Solve
[model,tlist,domainVolume]= fewcell.problemSetup(params,params.viz.showMesh);
fprintf('--> Solving...\n');
result= fewcell.solveProblem(model,tlist,params.solve);

%% Plot solution
% Prepare solution interpolation
fprintf('--> Interpolating solution...\n');
[AHLDistrib,x,y,interp_t,totalAHL]= ...
  fewcell.util.interpolateIntegrateAHL(model,result,params);
% Plot
fprintf('--> [paused] Press key to plot the solution\n'); pause;
fprintf('--> Plotting solution...\n');
fewcell.output(model,result,tlist,...
  AHLDistrib,x,y,interp_t,totalAHL,...
  params.viz.figID,params.viz.dynamicScaling);
