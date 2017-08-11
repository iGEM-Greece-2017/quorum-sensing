function [y,t]= run(y0,t,N)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
abstol= [1e-5,1e-6,1e-5,1e-3,1e-3,1e-3,1e-6,1e-6,1e-6];
opt= odeset('Jacobian', @singlecell.modelJacobian, ...
  'RelTol',1e-4,'AbsTol',abstol,'InitialStep',1);
[t,y]= ode23tb(@singlecell.model,t,y0,opt);
