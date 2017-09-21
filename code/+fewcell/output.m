function output(model,result,tlist, AHLDistrib,x,y,interp_t,totalAHL,params)

% Printed results
% [nM]=[nmol]/L => [mmol]= [nM]*[mm^3]
fprintf('Final total AHL at t=%.1fmin: %.3g [fmol]\n', interp_t(end), totalAHL(end));
fprintf('Final max [AHL] at t=%.1fmin: %.3g [nM]\n', interp_t(end), max(max(AHLDistrib(:,:,end))));
fprintf('Nodal max [AHL] at t=%.1fmin: %.3g [nM]\n', interp_t(end), max(result.NodalSolution(:,end)));

% Distribution plot
global distribAHL_interp_graphics;
if params.figID(1)>0
  distribAHL_interp_graphics= []; distribAHL_interp_graphics.first= true;
  for t= 1:length(interp_t)
    tic;
    fewcell.viz.plot(AHLDistrib,x,y,interp_t,totalAHL,params.figID(1),t,params.dynamicScaling);
    pauseTime= 1/22-toc;
    if pauseTime>0, pause(pauseTime); end
  end
end

% Singlecell plot
global yyResults;
if params.figID(2)>0
  singlecell.viz.plot(tlist, yyResults{1}, params.figID(2),params.logscaleSinglecell);
end

% Nodal solution
if params.figID(3)>0
  figure(params.figID(3)); clf;
  pdeplot(model,'XYData',result.NodalSolution(:,end), 'ZData',result.NodalSolution(:,end), 'Mesh','on','colormap','jet');
  view(0,90);
  nodalPlot= gca;
  nodalPlot.XLim= [0 params.domLim(1)];
  nodalPlot.YLim= [-params.domLim(2) 0];
end
