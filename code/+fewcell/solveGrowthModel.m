function [result,model]= solveGrowthModel(params,tlist,continued)
global growth;
global yyResults;
if nargin==3
  result= continued{1};
  model= continued{2};
  finalPointSaved= find(cellfun(@(x)~isempty(x),result(:,1)),1,'last');
  bactStart= (params.g.startRingIdx-1)*params.g.nLayers+1;
  yyResults= result{finalPointSaved,2}(bactStart: bactStart+length(yyResults)-1);
  bactEndState= cellfun(@(x) x(end,[1:5,7:9])', yyResults, 'UniformOutput',0); bactEndState= [bactEndState{:}];
else
  growth.tstep= [1;growth.tstep;length(tlist)];
  result= cell(length(growth.tstep)-1,2);
  model= cell(length(growth.tstep)-1,1);
  finalPointSaved= 1;
end
initConditions= 0;
fprintf('--> Solving...\n');
tic;
for i= finalPointSaved:length(growth.tstep)-1
  tlist_step= tlist(growth.tstep(i):growth.tstep(i+1));
  % Update both init conditions (solve.y0) and geometry (bactCenter0, nRings)
  if i>1
    params= fewcell.util.growthUpdateSolution(i,params,bactEndState);
    % function_handle to results_interpolant
    result_tend= util.sliceResult(model{i-1},result{i-1,1},length(result{i-1,1}.SolutionTimes));
    initConditions= @(loc) reshape(interpolateSolution(result_tend, loc.x,loc.y), size(loc.x));
  end
  % Create the new geometry
  tic;
  [params,model{i},domainVolume]= fewcell.problemSetup(params,initConditions);
  result{i,1}= fewcell.solveProblem(model{i},tlist_step,params);
  toc;
  bactEndState= cellfun(@(x) x(end,[1:5,7:9])', yyResults, 'UniformOutput',0); bactEndState= [bactEndState{:}];
  result{i,2}= cell(params.g.max_nRings*params.g.nLayers,1);
  bactStart= (params.g.startRingIdx-1)*params.g.nLayers+1;
  result{i,2}(1:bactStart-1)= {[zeros(length(tlist_step),10),ones(length(tlist_step),1)*length(yyResults)]};
  result{i,2}(bactStart: bactStart+length(yyResults)-1)= yyResults;
  result{i,2}(bactStart+length(yyResults):end)= {[zeros(length(tlist_step),10),ones(length(tlist_step),1)*length(yyResults)]};
  if nargin==2
    save(['data/tmpresults_',num2str(params.runID),'.mat'], 'params','model','result','tlist');
  else
    save(['data/tmpresults_',num2str(params.runID),'-continued.mat'], 'params','model','result','tlist');
  end
  fprintf('--> Growth step %d done\tt=%.0f/%dmin\n', i, tlist(growth.tstep(i+1)),tlist(end));
end
toc;
