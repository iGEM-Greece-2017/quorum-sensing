%% Runs the single cell model
t= linspace(0,60*25,400);
%{
global modelJacobian_minsvd;
modelJacobian_minsvd.svd= [];
modelJacobian_minsvd.stability= [];
modelJacobian_minsvd.j= [];
%}
first= true;
for N0= 10
  if ~first, pause; end
  fprintf('N0=%d\n',N0);
  y0= N0*[1.6;0;0;0;0;0;0;0;0;0;1];
  tic;
  y= singlecell.run_wMembrane(y0,t,N0);
  toc;

  if size(y,2)==10, y(:,10)= y(:,6); end
  singlecell.viz.plot(t,y,N0,1,false);
  singlecell.viz.plot(t,y,N0,2,true);
  fprintf('Nend: %.2g\n', y(end,11));
  fprintf('y1: %.2g\t y2: %.2g\t y6: %.2g\t y8: %.2g\t y9: %.2g\t y10: %.2g\n', ...
    y(end,1)/y(end,11),y(end,2)/y(end,11),y(end,6),y(end,8)/y(end,11),y(end,9)/y(end,11),y(end,10));
  %{
  figure(2); title([num2str(N0), ' cell']);
  yyaxis left; plot(modelJacobian_minsvd.svd(2:end));
  yyaxis right; plot(modelJacobian_minsvd.stability(2:end));
  legend('minsvd','svdratio');
  %}
  first= false;
end

save('code/+singlecell/+simplify/data', 'y','t');
