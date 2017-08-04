function [y,t]= run(y0,t,N)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
%abstol= [1e-18,1e-20,1e-19,1e-17,1e-18,1e-17,1e-25,1e-39,1e-49]*1e3;
abstol= 1e-6;
opt= odeset('Jacobian', @singlecell.modelJacobian, ...
  'RelTol',1e-5,'AbsTol',abstol,'InitialStep',1e-1);
[t,y]= ode23tb(@singlecell.model,t,y0,opt);
