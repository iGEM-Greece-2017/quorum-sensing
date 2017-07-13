%% Runs the single cell model
t= linspace(0,60000,10000);
y0= [1.6E-7;0;0;0;0;2e-5;0;0;0;0].*1e-9;
y= singlecell.run_wMembrane(y0,t);

singlecell.viz.plot(t,y0,y,1);
fprintf('y1: %.2g\t y2: %.2g\t y6: %.2g\t y8: %.2g\t y9: %.2g\n', ...
  y(end,1),y(end,2),y(end,6),y(end,8),y(end,9));
