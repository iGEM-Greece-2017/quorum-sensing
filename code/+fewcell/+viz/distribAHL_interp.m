function distribAHL_interp(u,x,y,total_u,times,domainLim,dynamicScale,figID,t)
global distribAHL_interp_graphics;

titleG1= {sprintf('t=%.2esec | step=%d',times(t),t), ...
  '[AHL] distribution (interpolated)'};
titleG2= {titleG1{1}, ...
  'Total [AHL]'};

if ~ishandle(figID) || distribAHL_interp_graphics.first
  distribAHL_interp_graphics.first= false;
  figure(figID); clf;
  fig= gcf; fig.KeyPressFcn= @util.graphicsPause;
  % [AHL] distribution
  subplot(1,3,[1 2]);
  imagesc(x,y, u(:,:,t));
  set(gca,'YDir','Normal');
  colormap hot;
  distribAHL_interp_graphics.c= colorbar;
  distribAHL_interp_graphics.c.Label.String= '[AHL]';
  if ~dynamicScale
    distribAHL_interp_graphics.globalScale= [min(u(:)) max(u(:))];
    %distribAHL_interp_graphics.c.Limits= distribAHL_interp_graphics.globalScale;
  end
  xlim([-1 1]*domainLim); ylim([-1 1]*domainLim*sin(pi/3));
  xlabel('x'); ylabel('y');
  title(titleG1);
  distribAHL_interp_graphics.g1= gca;
  
  % Total [AHL]
  subplot(1,3,3);
  plot(times(t),total_u(t), ...
    'LineStyle','-.','Marker','o','MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','r');
  xlim([times(1) times(end)]);
  ylim([min(total_u) max(total_u)]);
  grid minor; xlabel('time'); ylabel('[AHL]');
  title(titleG2);
  distribAHL_interp_graphics.g2= gca;
else
  % [AHL] distribution
  if dynamicScale
    distribAHL_interp_graphics.g1.CLimMode= 'auto';
    %distribAHL_interp_graphics.c.Limits= [min(u) max(u)];
  else
    %distribAHL_interp_graphics.c.Limits= distribAHL_interp_graphics.globalScale;
    distribAHL_interp_graphics.g1.CLim= distribAHL_interp_graphics.globalScale;
  end
  distribAHL_interp_graphics.g1.Children.CData= u(:,:,t);
  distribAHL_interp_graphics.g1.Title.String{1}= titleG1{1};
  
  % Total [AHL]
  distribAHL_interp_graphics.g2.Children.XData= ...
    [distribAHL_interp_graphics.g2.Children.XData, times(t)];
  distribAHL_interp_graphics.g2.Children.YData= ...
    [distribAHL_interp_graphics.g2.Children.YData, total_u(t)];
  distribAHL_interp_graphics.g2.Title.String{1}= titleG2{1};
end
drawnow;
end
