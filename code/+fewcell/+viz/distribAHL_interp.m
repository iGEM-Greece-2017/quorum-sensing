function distribAHL_interp(result,domainLim, resolution, dynamicScale, figID, t)
global distribAHL_interp_graphics;

titleG1= {sprintf('t=%.2esec | step=%d',result.SolutionTimes(t),t), ...
  '[AHL] distribution (interpolated)'};
titleG2= {titleG1{1}, ...
  'Total [AHL]'};

if ~ishandle(figID) || distribAHL_interp_graphics.first
  distribAHL_interp_graphics.first= false;
  
  % Precompute interpolation and integral for all times and store it
  x= linspace(-domainLim,domainLim, resolution);
  y= x*sin(pi/3); %y= flip(y);
  [X,Y]= meshgrid(x,y);
  u= @(x,y,t) reshape(interpolateSolution(result, x,y,t), size(x));
  distribAHL_interp_graphics.u= zeros([size(X),length(result.SolutionTimes)]);
  distribAHL_interp_graphics.total_u= zeros(1,length(result.SolutionTimes));
  for ti= t:length(result.SolutionTimes)
    distribAHL_interp_graphics.u(:,:,ti)= u(X,Y,ti);
    distribAHL_interp_graphics.total_u(ti)= integral2(@(x,y) u(x,y,ti),...
      -domainLim,domainLim, @(x)util.hexagonPerim(x,domainLim,-1), @(x)util.hexagonPerim(x,domainLim,1),...
      'AbsTol',0,'RelTol',4);
  end

  % Plot
  figure(figID); clf;
  fig= gcf; fig.KeyPressFcn= @util.graphicsPause;
  % [AHL] distribution
  subplot(1,3,[1 2]);
  imagesc(x,y, distribAHL_interp_graphics.u(:,:,t));
  set(gca,'YDir','Normal');
  colormap hot;
  distribAHL_interp_graphics.c= colorbar;
  distribAHL_interp_graphics.c.Label.String= '[AHL]';
  if ~dynamicScale
    distribAHL_interp_graphics.globalScale= [min(distribAHL_interp_graphics.u(:)) max(distribAHL_interp_graphics.u(:))];
    %distribAHL_interp_graphics.c.Limits= distribAHL_interp_graphics.globalScale;
  end
  xlim([-1 1]*domainLim); ylim([-1 1]*domainLim*sin(pi/3));
  xlabel('x'); ylabel('y');
  title(titleG1);
  distribAHL_interp_graphics.g1= gca;
  
  % Total [AHL]
  subplot(1,3,3);
  plot(result.SolutionTimes(t),distribAHL_interp_graphics.total_u(t), ...
    'LineStyle','-.','Marker','o','MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','r');
  xlim([result.SolutionTimes(1) result.SolutionTimes(end)]);
  ylim([min(distribAHL_interp_graphics.total_u) max(distribAHL_interp_graphics.total_u)]);
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
  distribAHL_interp_graphics.g1.Children.CData= distribAHL_interp_graphics.u(:,:,t);
  distribAHL_interp_graphics.g1.Title.String{1}= titleG1{1};
  
  % Total [AHL]
  distribAHL_interp_graphics.g2.Children.XData= ...
    [distribAHL_interp_graphics.g2.Children.XData, result.SolutionTimes(t)];
  distribAHL_interp_graphics.g2.Children.YData= ...
    [distribAHL_interp_graphics.g2.Children.YData, distribAHL_interp_graphics.total_u(t)];
  distribAHL_interp_graphics.g2.Title.String{1}= titleG2{1};
end
drawnow;
end
