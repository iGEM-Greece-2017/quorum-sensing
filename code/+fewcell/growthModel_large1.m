% Add growth model with evolving geometry to the basic distribModel
%
matlab_start;
clear params; clear global;
%% Parameters
% Units:
% Molarity: nM, Time: min, Length: mm, Quantity: fmol (nM*mm^3)
params.runID= randi(2147483644);
params.runName= "Small disk, growth, few growth.maxrings";
global enableSinglecellEq; enableSinglecellEq= true;  % false: debugging only
global enableGraphics; enableGraphics= true;
% time
params.t.tstop= 60*24;   % min
params.t.tstart= 0;
params.t.timePoints= 1200;
% coefficients
%params.c.c_agar= 7.1e-5*60;                 % [mm^2/min]
params.c.c_agar= 25e-5*60;                 % [mm^2/min]
params.c.c_cytoplasm= 35e-5*60;            % [mm^2/min]
params.c.d_AHL= 7e-5;                      % [1/min]
%params.c.d_AHL= 0.001;

% geometry
params.g.bactSize= 1e-3*[1,2.164];
params.g.init_bactCenter0= 1e-3*[120,-1.082];
params.g.max_nRings= 1652; params.g.nLayers= 2;
params.g.ringDist= 5;         % must be an odd number
params.g.layerSeparation= 1;
params.g.domainLim= [17,5.51];       % small disk
%params.g.domainLim= [1.7, .551];     % xs
%params.g.domainLim= [1.4, .45];     % xs
%params.g.domainLim= [.4,.1];         % tiny
%params.g.domainLim= [20,20]*1e-3;    % tiny

% Membrane r: 12.5mm -> (30oC)
% - Full nutrient agar: 10^10.1 cfu/membrane, growth.r=1.71
% *- 1/5 nutrient: 10^9.45 cfu/membrane, growth.r=1.5
% - 1/25 nutrient: 10^8.9 cfu/membrane, growth.r=1.5
% - no nutrients: 10^7.75 cfu/membrane, growth.r=.98

% growth
params.growth.on= true;    % enable/disable growth
% growth curve params
params.growth.gc.r= 1.5/60;
params.growth.gc.m= 0.52;
params.growth.gc.n= 3.5;
% growth step params
params.growth.r0= 2;        % How many rings of bacteria to start with
params.growth.dr= 6;        % How many rings of bacteria to add at each growth step
params.growth.min_dt= 35;
params.growth.maxRings= 300;

% mesh
params.m.Hgrad= 1.5;
params.m.HmaxCoeff= 1/12;
% init
params.solve.y0= [1.5347;0;0;0;0; 0;0;0];  % [nM]
% solve
params.solve.AbsTol_y= [1e-3,2e-3,2e-3,5e-3,5e-3,  2e-3,2e-3,1e-3]*1e0;
params.solve.AbsTol= 1e-4;    % for diffusion nodes
params.solve.RelTol= 1e-3;
params.solve.FeatureSize= min(params.g.bactSize);
params.solve.reportStatistics= 'on';

% viz
params.viz.showMesh= true;
params.viz.showGrowthCurve= true;
params.viz.domLim= [[0;0],params.g.domainLim'./1];
params.viz.interpResolution= 120;
params.viz.integrateAbstol= 1;
params.viz.timePoints= floor(params.t.timePoints/6);
params.viz.dynamicScaling= true;
params.viz.logscaleSinglecell= false;
params.viz.figID= [5,6,0];

%% Solve
[params,tlist,bactRingDensity]= fewcell.initSetup(params);
[result,model]= fewcell.solveGrowthModel(params,tlist);
finalPointSaved= find(cellfun(@(x)~isempty(x),result(:,1)),1,'last');
%}
%save(['data/tmpresults_',num2str(params.runID),'.mat'], 'params','model','result','tlist');
%% Plot solution
% Prepare solution interpolatms.g.ringDist= 1*1e-3; params.g.layerSeparation= 1.082*1e-3;
fprintf('--> Interpolating solution...\n');
AHLDistrib= zeros(params.viz.interpResolution);
interp_t= []; totalAHL= [];
totalVizTimepoints= params.viz.timePoints;
global growth;
global yyResults;
yyResults= cellfun(@(x) x(1,:), result{1,2}, 'uniformoutput',0);
tic;
for i= 1:finalPointSaved
  tstep= 1:growth.tstep(i+1)-growth.tstep(i);
  tpart= (tlist(growth.tstep(i+1))-tlist(growth.tstep(i)))/(tlist(end)-tlist(1));
  params.viz.timePoints= ceil(totalVizTimepoints*tpart);
  [tAHLDistrib,x,y,tInterp_t,ttotalAHL]= ...
    fewcell.util.interpolateIntegrateAHL(model{i},result{i,1},params,tstep(1),tstep(end));
  params.viz.timePoints= totalVizTimepoints;
  AHLDistrib(:,:,size(AHLDistrib,3)+1:size(tAHLDistrib,3)+size(AHLDistrib,3))= tAHLDistrib;
  interp_t= [interp_t,tInterp_t];
  totalAHL= [totalAHL,ttotalAHL];
  yyResults= cellfun(@(x,y) [x;y(2:end,:)], yyResults,result{i,2}, 'uniformoutput',0);
end
AHLDistrib= AHLDistrib(:,:,2:end);
toc;
% Plot
fprintf('--> Plotting solution...\n');
fewcell.output(model{i},result{i,1},tlist(growth.tstep(1):growth.tstep(i+1)),AHLDistrib,x,y,interp_t,totalAHL,params.viz);
