%% Runs the single cell model
t= linspace(0,60*50,900);
first= true;
for N0= 5e4
  if ~first, pause; end
  fprintf('N0=%d\n',N0);
  y0= [1.5347;0;0;0;0;0;0;0;0;0;N0];
  tic;
  y= singlecell.run_weber(y0,t,true,true);
  assert(length(y)==length(t));
  toc;

  if size(y,2)==10, y(:,10)= y(:,6); end
  singlecell.viz.plot(t,y,1,true);
  fprintf('Final Ncell: %.2g\n', y(end,11));
  fprintf('y1: %.2g\t y2: %.2g\t y4: %.2g\t y6: %.2g\t y8: %.2g\t y9: %.2g\t y10: %.2g\n', ...
    y(end,1),y(end,2),y(end,4),y(end,6),y(end,8),y(end,9),y(end,10));
  first= false;
end

save('code/+singlecell/+simplify/data', 'y','t');
