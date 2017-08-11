%% Runs the single cell model
t= linspace(0,60*80,400);
%t= linspace(0,300,300);
%spaceVol= 1^3; %mm3
spaceVol= 1.26^3*1e-9; % no dilution
%{
global modelJacobian_minsvd;
modelJacobian_minsvd.svd= [];
modelJacobian_minsvd.stability= [];
modelJacobian_minsvd.j= [];
%}
first= true;
for N= logspace(0,6,4)
  if ~first, pause; end
  fprintf('N=%d\n',N);
  y0= N*[1.6;0;0;0;0;0;0;0;0;0] .* singlecell.yNormalization(1:10)';
  tic;
  y= singlecell.run_wMembrane(y0,t,N,spaceVol);
  toc;
  y= y ./ repmat(singlecell.yNormalization(1:10),size(y,1),1);

  if size(y,2)==9, y(:,10)= y(:,6); end
  singlecell.viz.plot(t,y,N,1);
  fprintf('y1: %.2g\t y2: %.2g\t y6: %.2g\t y8: %.2g\t y9: %.2g\t y10: %.2g\n', ...
    y(end,1)/N,y(end,2)/N,y(end,6),y(end,8)/N,y(end,9)/N,y(end,10));
  %{
  figure(2); title([num2str(N), ' cell']);
  yyaxis left; plot(modelJacobian_minsvd.svd(2:end));
  yyaxis right; plot(modelJacobian_minsvd.stability(2:end));
  legend('minsvd','svdratio');
  %}
  first= false;
end

save('code/+singlecell/+simplify/data', 'y','t');
