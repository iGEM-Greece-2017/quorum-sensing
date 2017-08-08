% Estimate simplified model in full singlecell
clear;
%% Some parameters, data
N= 100;
x0= N*[1.6;0;0;0;0;0;0;0;0;0];
%abstol= [1e-5,1e-6,1e-5,1e-3,1e-3,1e-3,1e-6,1e-6,1e-6,1e-5];
nt= 20;
ts= 4.5;

k1=0.01; k2=1; k3=0.05; k4=1;  dS=0.002; dSS=0.002;
load('code/+singlecell/+simplify/estimatedModel');  % newModel
paramVal= newModel.Report.Parameters.ParVector;
%% Model setup
% Original model
singlecellModel= idnlgrey(@singlecell.simplify.origModel_full,[2 0 10], [],x0, 0,'TimeUnit','minutes');
singlecellModel.FileArgument= {N};
singlecellModel.SimulationOptions.Solver= 'ode15s';
singlecellModel.SimulationOptions.AbsTol= 1e-3;
singlecellModel.SimulationOptions.RelTol= 1e-3;

% New model
params= struct('Name', {'k1','k2','a','order_f','order_b'},...
  'Value', {paramVal(1), paramVal(2), paramVal(3), paramVal(4), paramVal(5)},...
  'Minimum', {0,0,-dSS,0.5,0.5},...
  'Maximum', {10*max([k1,k3,k1*k3]), 10*max([k2,k4,k2*k4]), 5*dS, 1.5,1.5},...
  'Unit', '', 'Fixed', {0,0,0,0,0});
newModel= idnlgrey(@singlecell.simplify.newModel_full,[2 0 9], params,...
  [x0(1:6);x0(8:10)], 0,'TimeUnit','minutes');
newModel.FileArgument= {N};
newModel.SimulationOptions.Solver= 'ode15s';
newModel.SimulationOptions.AbsTol= 1e-3;
newModel.SimulationOptions.RelTol= 1e-3;

% Simulate original
zData= iddata(zeros(nt,2),[], ts, 'TimeUnit','minutes');
singlecellData= sim(singlecellModel,zData);
newInitData= sim(newModel,zData);
singlecellData.Ts= ts;
% Compare initial newModel with singlecell
figure(3); clf;
%
subplot(211); compare(singlecellData,newModel); grid minor; title('Initial-Singlecell');
%}

%% Estimate (on "orig" data)
opt= nlgreyestOptions('Display','on');
opt.Regularization.Nominal= 'model';
opt.Regularization.Lambda= 2;
%opt.Regularization.R= [1,1,1,100];
opt.SearchOption.MaxIter= 50;
newModel= nlgreyest(singlecellData,newModel,opt);
% Compare newModel with origModel
subplot(212); compare(singlecellData,newModel); grid minor;
title('Estimated-Singlecell');

%% Save estimated model
newModel.Report.Parameters.ParVector
save('code/+singlecell/+simplify/estimatedModel_full', 'newModel');
