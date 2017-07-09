function [y,t]= run_wMembrane(y0,t)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
abstol= [1e-18,1e-20,1e-19,1e-17,1e-18,1e-17,1e-25,1e-39,1e-49,1e-17];
opt= odeset('Jacobian', @singlecell.modelJacobian_wMembrane, ...
  'RelTol',1e-5,'AbsTol',abstol,'InitialStep',1e-5);
[t,y]= ode23tb(@singlecell.model_wMembrane,t,y0,opt);
