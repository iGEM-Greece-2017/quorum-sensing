function distribAHL_node(model,result,umin,domainLim, dynamicScale, figID, t)
global distribAHL_node_graphics;

u= result.NodalSolution(:,t); u= u+umin;  % Can't have negative concentrations, so shift all values slightly up
u= log10(u);
%u= log10(u);
if ~ishandle(figID) || t==1
  figure(figID); clf;
  fig= gcf; fig.KeyPressFcn= @util.graphicsPause;
  
  % [AHL] distribution
  subplot(1,3,[1 2]);
  pdeplot(model,'XYData',u, ...
    'ColorMap','hot','Mesh','on');
  distribAHL_node_graphics.c= colorbar;
  distribAHL_node_graphics.c.Label.String= 'log10([AHL])';
  if ~dynamicScale
    uAll= result.NodalSolution(:); uAll= uAll+umin;
    uAll= log10(uAll);
    distribAHL_node_graphics.globalScale= [min(uAll) max(uAll)];
    distribAHL_node_graphics.c.Limits= distribAHL_node_graphics.globalScale;
  end
  xlim([-1 1]*domainLim); ylim([-1 1]*domainLim);
  title({['t=',num2str(result.SolutionTimes(t)), ' step=',num2str(t)], '[AHL] distribution'});
  distribAHL_node_graphics.g1= gca;
  
  % Total [AHL]
  subplot(1,3,3);
  plot(result.SolutionTimes(t),sum(10.^u), ...
    'LineStyle','-.','Marker','o','MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','r');
  xlim([result.SolutionTimes(1) result.SolutionTimes(end)]);
  grid minor; xlabel('time'); ylabel('[AHL]');
  title({['t=',num2str(result.SolutionTimes(t)), ' step=',num2str(t)], 'Total [AHL] at nodes'});
  distribAHL_node_graphics.g2= gca;
else
  if dynamicScale
    distribAHL_node_graphics.g1.CLimMode= 'auto';
    %distribAHL_node_graphics.c.Limits= [min(u) max(u)];
  else
    %distribAHL_node_graphics.c.Limits= distribAHL_node_graphics.globalScale;
    distribAHL_node_graphics.g1.CLim= distribAHL_node_graphics.globalScale;
  end
  distribAHL_node_graphics.g1.Children.CData= u(result.Mesh.Elements);
  distribAHL_node_graphics.g1.Title.String= ['t=',num2str(result.SolutionTimes(t)), ' step=',num2str(t)];
  
  distribAHL_node_graphics.g2.Children.XData= [distribAHL_node_graphics.g2.Children.XData, result.SolutionTimes(t)];
  distribAHL_node_graphics.g2.Children.YData= [distribAHL_node_graphics.g2.Children.YData, sum(10.^u)];
end
drawnow;
end
