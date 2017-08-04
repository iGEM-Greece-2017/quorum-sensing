function output(AHLDistrib,x,y,interp_t,totalAHL,figID,dynamicScaling)

% Printed results
fprintf('Final total AHL at t=%.1fsec: %.3g [mol]\n', interp_t(end), totalAHL(end));
fprintf('Final max [AHL] at t=%.1fsec: %.3g [M]\n\n', interp_t(end), max(max(AHLDistrib(:,:,end))));
global distribAHL_interp_graphics;
distribAHL_interp_graphics= []; distribAHL_interp_graphics.first= true;
for t= 1:length(interp_t)
  fewcell.viz.plot(AHLDistrib,x,y,interp_t,totalAHL,figID(1),t,dynamicScaling);
end
end
