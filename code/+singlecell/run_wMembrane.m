function [y,t]= run_wMembrane(y0,t,N,spaceVol)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
abstol= [1e-5,1e-6,1e-5,1e-3,1e-3,1e-3,1e-6,1e-6,1e-6,1e-5];
opt= odeset('Jacobian', @(t,y)singlecell.modelJacobian_wMembrane(t,y,N,spaceVol), ...
  'RelTol',1e-4,'AbsTol',abstol,'InitialStep',1);
[t,y]= ode15s(@(t,y)singlecell.model_wMembrane(t,y,N,spaceVol),t,y0,opt);
