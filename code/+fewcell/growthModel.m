% Add growth model with evolving geometry to the basic distribModel
matlab_start;
clear params; clear global;
%% Parameters
% Units:
% Molarity: nM, Time: min, Length: mm, Quantity: fmol (nM*mm^3)
params.runID= randi(2147483644);
global enableSinglecellEq; enableSinglecellEq= true;  % false: debugging only
global enableGraphics; enableGraphics= true;
% time
params.t.tstop= 60*18;   % min
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
params.g.init_bactCenter0= 1e-3*[300,-1.082];
params.g.max_nRings= 50; params.g.nLayers= 2;
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
params.growth.on= true;    % enable/disable growth
% growth curve params
params.growth.gc.r= 1.5/60;
params.growth.gc.m= 0.52;
params.growth.gc.n= 3.5;
% growth step params
params.growth.r0= 2;      % How many rings of bacteria to start with
params.growth.dr= 8;      % How many rings of bacteria to add at each growth step
params.growth.maxRings= 80;

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
params.viz.figID= [5,6,0];

%% Solve
global growth;
global yyResults;
[params,tlist,bactRingDensity]= fewcell.initSetup(params);
growth.tstep= [1;growth.tstep];
initConditions= 0;
result= cell(length(growth.tstep)-1,3);
model= cell(length(growth.tstep)-1,1);
fprintf('--> Solving...\n');
tic;
for i= 1:length(growth.tstep)-1
  tlist_step= tlist(growth.tstep(i):growth.tstep(i+1));
  % Update both init conditions (solve.y0) and geometry (bactCenter0, nRings)
  if i>1
    params= fewcell.util.growthUpdateSolution(i,params,result(i-1,:));
    % function_handle to results_interpolant
    result_tend= util.sliceResult(model{i-1},result{i-1,1},length(result{i-1,1}.SolutionTimes));
    initConditions= @(loc) reshape(interpolateSolution(result_tend, loc.x,loc.y), size(loc.x));
  end
  % Create the new geometry
  [params,model{i},domainVolume]= fewcell.problemSetup(params,initConditions);
  result{i,1}= fewcell.solveProblem(model{i},tlist_step,params);
  result{i,2}= cellfun(@(x) x(end,[1:5,7:9])', yyResults, 'UniformOutput',0); result{i,2}= [result{i,2}{:}];
  result{i,3}= cell(params.g.max_nRings*params.g.nLayers,1);
  result{i,3}(1:length(yyResults))= yyResults;
  result{i,3}(length(yyResults)+1:end)= {[zeros(length(tlist_step),10),ones(length(tlist_step),1)*length(yyResults)]};
  fprintf('--> Growth step %d done\tt=%.0f/%dsec\n', i, tlist(growth.tstep(i)),tlist(end));
end

toc;
save(['data/tmpresults_',num2str(params.runID),'.mat'], 'params','model','result','tlist');
%% Plot solution
% Prepare solution interpolatms.g.ringDist= 1*1e-3; params.g.layerSeparation= 1.082*1e-3;
fprintf('--> Interpolating solution...\n');
AHLDistrib= zeros(params.viz.interpResolution);
interp_t= []; totalAHL= [];
totalVizTimepoints= params.viz.timePoints;
tic;
for i= 1:length(growth.tstep)-1
  tstep= 1:growth.tstep(i+1)-growth.tstep(i);
  tpart= (tlist(growth.tstep(i+1))-tlist(growth.tstep(i)))/(tlist(end)-tlist(1));
  params.viz.timePoints= ceil(totalVizTimepoints*tpart);
  [tAHLDistrib,x,y,tInterp_t,ttotalAHL]= ...
    fewcell.util.interpolateIntegrateAHL(model{i},result{i,1},params,tstep(1),tstep(end));
  params.viz.timePoints= totalVizTimepoints;
  AHLDistrib(:,:,size(AHLDistrib,3)+1:size(tAHLDistrib,3)+size(AHLDistrib,3))= tAHLDistrib;
  interp_t= [interp_t,tInterp_t];
  totalAHL= [totalAHL,ttotalAHL];
end
AHLDistrib= AHLDistrib(:,:,2:end);
toc;
% Plot
fprintf('--> Plotting solution...\n');
fewcell.output(model{end},result{end,1},tlist(growth.tstep(end-1):growth.tstep(end)),AHLDistrib,x,y,interp_t,totalAHL,params.viz);
