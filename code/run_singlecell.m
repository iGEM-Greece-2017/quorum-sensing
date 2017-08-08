%% Runs the single cell model
t= linspace(0,60*30,400);
%t= linspace(0,300,300);
%{
y0= [1.6E-7;0;0;0;0;2e-5;0;0;0;0].*1e-9;
y= singlecell.run_wMembrane(y0,t);
%}
global modelJacobian_minsvd;
modelJacobian_minsvd.svd= [];
modelJacobian_minsvd.stability= [];
modelJacobian_minsvd.j= [];
first= true;
for N= 5%logspace(0,6,4)
  if ~first, pause; end
  fprintf('N=%d\n',N);
  y0= N*[1.6;0;0;0;0;0;0;0;0;0] .* singlecell.yNormalization(1:10)';
  tic;
  y= singlecell.run_wMembrane(y0,t,N);
  toc;
  y= y ./ repmat(singlecell.yNormalization(1:10),size(y,1),1);

  singlecell.viz.plot(t,y,N,1);

  if size(y,2)==9, y(:,10)= y(:,6); end
  fprintf('y1: %.2g\t y2: %.2g\t y6: %.2g\t y8: %.2g\t y9: %.2g\t y10: %.2g\n', ...
    y(end,1),y(end,2),y(end,6),y(end,8),y(end,9),y(end,10));
  %{
  figure(2); title([num2str(N), ' cell']);
  yyaxis left; plot(modelJacobian_minsvd.svd(2:end));
  yyaxis right; plot(modelJacobian_minsvd.stability(2:end));
  legend('minsvd','svdratio');
  %}
  first= false;
end

save('code/+singlecell/+simplify/data', 'y','t');
