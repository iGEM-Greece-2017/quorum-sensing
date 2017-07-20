function distribAHL_interp(u,x,y,total_u,times,zoomin,dynamicScale,figID,t)
global distribAHL_interp_graphics;

titleG1= {sprintf('t=%.2esec | step=%d',times(t),t), ...
  '[AHL] distribution'};
titleG2= {titleG1{1}, ...
  'Total AHL'};
titleG3= {titleG1{1}, ...
  'Max [AHL]'};

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
  distribAHL_interp_graphics.c.Label.String= '[AHL] [M]';
  if ~dynamicScale
    distribAHL_interp_graphics.globalScale= [min(u(:,:,t:end)) max(u(:,:,t:end))+eps(max(u(:,:,t:end)))];
  end
  axis tight;
  xlabel('x [m]'); ylabel('y [m]');
  title(titleG1);
  distribAHL_interp_graphics.g1= gca;
  
  % Total AHL
  subplot(2,3,3);
  plot(times(t),total_u(t), ...
    'LineStyle','-.','Marker','o','MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','r');
  xlim([times(t) times(end)]);
  ylim([min(total_u(t:end)) max(total_u)+eps(max(total_u(t:end)))]);
  grid minor; xlabel('time [s]'); ylabel('AHL [mol]');
  title(titleG2);
  distribAHL_interp_graphics.g2= gca;
  
  % Max [AHL]
  subplot(2,3,6);
  plot(times(t),max(max(u(:,:,t))), ...
    'LineStyle','-.','Marker','d','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','r');
  xlim([times(t) times(end)]);
  ylim([min(max(max(u(:,:,t:end),[],1),[],2)), ...
        max(max(max(u(:,:,t:end),[],1),[],2)) + ...
        eps(max(max(max(u(:,:,t:end),[],1),[],2)))]);
  grid minor; xlabel('time [s]'); ylabel('[AHL] [M]');
  title(titleG3);
  distribAHL_interp_graphics.g3= gca;
else
  % [AHL] distribution
  ut= u(:,:,t);
  distribAHL_interp_graphics.g1.Children.CData= ut;
  distribAHL_interp_graphics.g1.Title.String{1}= titleG1{1};
  if dynamicScale
    distribAHL_interp_graphics.g1.CLim= [min(ut(:)) max(ut(:))];
  else
    %distribAHL_interp_graphics.c.Limits= distribAHL_interp_graphics.globalScale;
    distribAHL_interp_graphics.g1.CLim= distribAHL_interp_graphics.globalScale;
  end
  
  % Total AHL
  distribAHL_interp_graphics.g2.Children.XData= ...
    [distribAHL_interp_graphics.g2.Children.XData, times(t)];
  distribAHL_interp_graphics.g2.Children.YData= ...
    [distribAHL_interp_graphics.g2.Children.YData, total_u(t)];
  distribAHL_interp_graphics.g2.Title.String{1}= titleG2{1};
  
  % Max [AHL]
  distribAHL_interp_graphics.g3.Children.XData= ...
    [distribAHL_interp_graphics.g3.Children.XData, times(t)];
  distribAHL_interp_graphics.g3.Children.YData= ...
    [distribAHL_interp_graphics.g3.Children.YData, max(max(u(:,:,t)))];
  distribAHL_interp_graphics.g3.Title.String{1}= titleG3{1};
end
drawnow;
end
