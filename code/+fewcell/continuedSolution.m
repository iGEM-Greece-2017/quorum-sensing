% Continue solving a previous solution
matlab_start;
%clear params; clear global;
%% Parameters
%load('data/results/final/tmpresults_677792149.mat');
%load('data/tmpresults_353316463.mat');
load('data/tmpresults_1527163184.mat');
global enableSinglecellEq; enableSinglecellEq= true;  % false: debugging only
global enableGraphics; enableGraphics= true;

params.solve.AbsTol_y= [1e-3,1e-2,1e-2,1e-3,1e-3,  1e-6,1e-6,1e-4]*1e-1;
params.solve.AbsTol= 1e-8;    % for diffusion nodes
params.solve.RelTol= 1e-2;

%% Solve
global growth;
global yyResults;
try growth= growthLoc; end
[result,model]= fewcell.solveGrowthModel(params,tlist,{result,model});
finalPointSaved= find(cellfun(@(x)~isempty(x),result(:,1)),1,'last');
%}
%save(['data/tmpresults_',num2str(params.runID),'.mat'], 'params','model','result','tlist');
%% Plot solution
% Prepare solution interpolatms.g.ringDist= 1*1e-3; params.g.layerSeparation= 1.082*1e-3;
fprintf('--> Interpolating solution...\n');
AHLDistrib= zeros(params.viz.interpResolution);
interp_t= []; totalAHL= [];
totalVizTimepoints= params.viz.timePoints;
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
