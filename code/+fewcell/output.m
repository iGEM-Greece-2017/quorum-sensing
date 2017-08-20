function output(model,result,tlist, AHLDistrib,x,y,interp_t,totalAHL,figID,dynamicScaling)

% Printed results
% [nM]=[nmol]/L => [mmol]= [nM]*[mm^3]
fprintf('Final total AHL at t=%.1fmin: %.3g [fmol]\n', interp_t(end), totalAHL(end));
fprintf('Final max [AHL] at t=%.1fmin: %.3g [nM]\n', interp_t(end), max(max(AHLDistrib(:,:,end))));
fprintf('Nodal max [AHL] at t=%.1fmin: %.3g [nM]\n', interp_t(end), max(result.NodalSolution(:,end)));

% Distribution plot
global distribAHL_interp_graphics;
if figID(1)>0
  distribAHL_interp_graphics= []; distribAHL_interp_graphics.first= true;
  for t= 1:length(interp_t)
    fewcell.viz.plot(AHLDistrib,x,y,interp_t,totalAHL,figID(1),t,dynamicScaling);
  end
end

% Singlecell plot
global yyResults;
if figID(2)>0
  singlecell.viz.plot(tlist, yyResults{1}, 1,figID(2));
end

% Nodal solution
if figID(3)>0
  figure(figID(3));
  pdeplot(model,'XYData',result.NodalSolution(:,end),'Mesh','on','colormap','hot');
  nodalPlot= gca;
  nodalPlot.XLim= [0 0.03];
  nodalPlot.YLim= [-0.015 0];
end
