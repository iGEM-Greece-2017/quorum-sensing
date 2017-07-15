function [y,t]= run_wMembrane(y0,t,N)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
abstol= [1e-18,1e-20,1e-19,1e-17,1e-18,1e-17,1e-25,1e-39,1e-49,1e-17];
opt= odeset('Jacobian', @(t,y)singlecell.modelJacobian_wMembrane(t,y,N), ...
  'RelTol',1e-8,'AbsTol',abstol,'InitialStep',1e-5);
[t,y]= ode15s(@(t,y)singlecell.model_wMembrane(t,y,N),t,y0,opt);
