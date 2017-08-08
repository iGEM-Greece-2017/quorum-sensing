% Estimate simplified model
clear;
%% Some parameters, data
% The kinetic constants are in (1/min, 1/(min*nM))
k1=0.01; k2=1; k3=0.05; k4=1;  dS=0.002; dSS=0.002;
load('code/+singlecell/+simplify/data');  % y,t

%% Model setup
% Original model
origModel= idnlgrey(@singlecell.simplify.origModel_part,[2 4 2], [],[0;0], 0,'TimeUnit','minutes');
origModel.SimulationOptions.Solver= 'ode15s';

% New model
params= struct('Name', {'k1','k2','a','order'},...
  'Value', {min(k1,k3), min(k2,k4), 2*dS, 1.03},...
  'Minimum', {0,0,0,0.5},...
  'Maximum', {10*max([k1,k3,k1*k3]), 10*max([k2,k4,k2*k4]), 5*dS, 1.5},...
  'Unit', '', 'Fixed', {0,0,0,0});
newModel= idnlgrey(@singlecell.simplify.newModel_part,[2 4 1], params,0, 0,'TimeUnit','minutes');

% Simulate original
y= y(:,:);
data= iddata([y(:,8), -k1*y(:,5).*y(:,6)+k2*y(:,7)], ...
             y(:,[5,6,1,9]), t(2)-t(1), ...
             'TimeUnit','minutes', 'InterSample',{'foh';'foh';'foh';'foh'});
origData= sim(origModel,data);
origData.InputData= y(:,[5,6,1,9]);
origData.InterSample= {'foh';'foh';'foh';'foh'}; origData.Ts= t(2)-t(1);
% Compare initial newModel with origModel
figure(2); clf;
%
%subplot(221); compare(data,origModel); grid minor; title('Original-Singlecell');
%subplot(223); compare(origData,newModel); grid minor; title('Initial-Original');
%}

%% Estimate (on "orig" data)
opt= nlgreyestOptions('Display','on');
opt.Regularization.Nominal= 'model';
opt.Regularization.Lambda= 1;
opt.Regularization.R= [50,1,1,100];
opt.SearchOption.MaxIter= 50;
opt.SearchOption.TolFun= 1e-7;
opt.SearchOption.Advanced.UseParallel= 1;
newModel= nlgreyest(origData,newModel,opt);
% Compare newModel with origModel/singlecell data
subplot(222); compare(data,newModel); grid minor; title('Estimated-Singlecell');
%subplot(224); compare(origData,newModel); grid minor; title('Estimated-Original');

%% Save estimated model
newModel.Report.Parameters.ParVector
save('code/+singlecell/+simplify/estimatedModel', 'newModel');
