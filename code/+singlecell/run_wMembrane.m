function [y,t]= run_wMembrane(y0,t,N)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
abstol= [1e-12,1e-14,1e-13,1e-10,1e-10,1e-15,1e-15,1e-19,1e-21,1e-15]*1e-2;
opt= odeset('Jacobian', @(t,y)singlecell.modelJacobian_wMembrane(t,y,N), ...
  'RelTol',1e-7,'AbsTol',abstol,'InitialStep',1e-2);
[t,y]= ode15s(@(t,y)singlecell.model_wMembrane(t,y,N),t,y0,opt);
