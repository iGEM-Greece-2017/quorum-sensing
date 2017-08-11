function plot(u,x,y,times,total_u,figID,t,dynamicScale)

times_h= times/60;
global distribAHL_interp_graphics;
if nargin<8, dynamicScale= true; end
titles{1}= {sprintf('t=%.2ehours | step=%d',times_h(t),t), '[AHL] distribution'};
titles{2}= {'Total AHL (entire domain)'};
titles{3}= {'Max [AHL] (zoomed domain)'};
% On the first call, create the graphics
if ~ishandle(figID) || distribAHL_interp_graphics.first
  createGraphics(u,x,y,times_h,total_u,figID,t,titles,dynamicScale);
else  % On subsequent calls, only update them
  updateGraphics(u,total_u,times_h,t,titles,dynamicScale);
end
drawnow;
end

function createGraphics(u,x,y,times_h,total_u,figID,t,titles,dynamicScale)
  global distribAHL_interp_graphics;
  distribAHL_interp_graphics.first= false;
  figure(figID); clf;
  fig= gcf; fig.KeyPressFcn= @util.graphicsPause;
  % [AHL] distribution
  subplot(1,3,[1 2]);
  imagesc(x,y, u(:,:,t));
  set(gca,'YDir','Normal');
  colormap hot;
  distribAHL_interp_graphics.c= colorbar;
  distribAHL_interp_graphics.c.Label.String= '[AHL] [nM]';
  if ~dynamicScale
    ut= u(:,:,t:end);
    distribAHL_interp_graphics.globalScale= [min(ut(:)) max(ut(:))+eps(max(ut(:)))];
  end
  axis tight;
  xlabel('x [mm]'); ylabel('y [mm]');
  title(titles{1});
  distribAHL_interp_graphics.g1= gca;
  
  % Total AHL
  subplot(2,3,3);
  plot(times_h(t),total_u(t), ...
    'LineStyle','-.','Marker','o','MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','r');
  xlim([times_h(t) times_h(end)]);
  ylim([min(total_u(t:end)) max(total_u)+eps(max(total_u(t:end)))]);
  grid minor; xlabel('time [hour]'); ylabel('AHL [fmol]');
  title(titles{2});
  distribAHL_interp_graphics.g2= gca;
  
  % Max [AHL]
  subplot(2,3,6);
  plot(times_h(t),max(max(u(:,:,t))), ...
    'LineStyle','-.','Marker','d','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','r');
  xlim([times_h(t) times_h(end)]);
  ylim([min(max(max(u(:,:,t:end),[],1),[],2)), ...
        max(max(max(u(:,:,t:end),[],1),[],2)) + ...
        eps(max(max(max(u(:,:,t:end),[],1),[],2)))]);
  grid minor; xlabel('time [hour]'); ylabel('[AHL] [nM]');
  title(titles{3});
  distribAHL_interp_graphics.g3= gca;
end

function updateGraphics(u,total_u,times_h,t,titles,dynamicScale)
  global distribAHL_interp_graphics;
  % [AHL] distribution
  ut= u(:,:,t);
  distribAHL_interp_graphics.g1.Children.CData= ut;
  distribAHL_interp_graphics.g1.Title.String{1}= titles{1}{1};
  if dynamicScale
    distribAHL_interp_graphics.g1.CLim= [min(ut(:)) max(ut(:))];
  else
    %distribAHL_interp_graphics.c.Limits= distribAHL_interp_graphics.globalScale;
    distribAHL_interp_graphics.g1.CLim= distribAHL_interp_graphics.globalScale;
  end
  
  % Total AHL
  distribAHL_interp_graphics.g2.Children.XData= ...
    [distribAHL_interp_graphics.g2.Children.XData, times_h(t)];
  distribAHL_interp_graphics.g2.Children.YData= ...
    [distribAHL_interp_graphics.g2.Children.YData, total_u(t)];
  
  % Max [AHL]
  distribAHL_interp_graphics.g3.Children.XData= ...
    [distribAHL_interp_graphics.g3.Children.XData, times_h(t)];
  distribAHL_interp_graphics.g3.Children.YData= ...
    [distribAHL_interp_graphics.g3.Children.YData, max(max(u(:,:,t)))];
end
