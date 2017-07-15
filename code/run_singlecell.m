%% Runs the single cell model
t= linspace(0,90000,5000);
%{
y0= [1.6E-7;0;0;0;0;2e-5;0;0;0;0].*1e-9;
y= singlecell.run_wMembrane(y0,t);
%}
for N= [1 1e2 1e4 1e6 1e7 1e8]
fprintf('N=%d\n',N);
y0= N*[1.6E-7;0;0;0;0;0;0;0;0;0].*1e-9;
y= singlecell.run_wMembrane(y0,t,N);

singlecell.viz.plot(t,y0,y,N,1);
fprintf('y1: %.2g\t y2: %.2g\t y6: %.2g\t y8: %.2g\t y9: %.2g\t y10: %.2g\n', ...
  y(end,1),y(end,2),y(end,6),y(end,8),y(end,9),y(end,10));
pause;
end
