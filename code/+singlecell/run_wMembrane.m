function [y,t]= run_wMembrane(y0,t,N)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
%abstol= [1e-10,5e-12,1e-11,1e-9,1e-8, 1e-13, 1e-13,1e-16,1e-18 ,1e-15]*1e0;
%abstol= [1e-4,1e-2,1e-2,1e-9,1e-8, 1e-13, 1e-13,1e-16,1e-18 ,1e-15]*1e0;
abstol= 1e-5;
opt= odeset('Jacobian', @(t,y)singlecell.modelJacobian_wMembrane(t,y,N), ...
  'RelTol',1e-5,'AbsTol',abstol,'InitialStep',1e-1);
[t,y]= ode15s(@(t,y)singlecell.model_wMembrane(t,y,N),t,y0,opt);
