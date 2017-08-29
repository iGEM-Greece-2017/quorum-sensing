function [y,t]= run_weber(y0,t,growthOn)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
abstol= [1e-5,1e-6,1e-5,1e-3,1e-3,1e-3,1e-6,1e-6,1e-6,1e-5,1e-3];
opt= odeset('Jacobian', [],...
  'RelTol',1e-4,'AbsTol',abstol,'InitialStep',1e-3);
[t,y]= ode15s(@(t,y)singlecell.model_weber(t,y,growthOn,true),t,y0,opt);
