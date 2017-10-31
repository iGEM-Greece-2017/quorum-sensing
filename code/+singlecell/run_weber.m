function [y,t]= run_weber(y0,t,growthOn,diffusionOn)
% Run the singlecell model. Given a starting state and start/end times, returns the model's
% state and their corresponding times
%abstol= [1e-5,1e-6,1e-5,1e-4,1e-4,1e-3,1e-6,1e-6,1e-6,1e-5,1e-3]*1e0;
abstol= [1e-3,2e-3,2e-3,5e-3,5e-3, 2e-3, 2e-3,2e-3,1e-3, 2e-3,1e-1]*1e0;

  %abstol= [];
%abstol= 1e-12;
opt= odeset('Jacobian', @(t,y)singlecell.modelJacobian_weberWrapper(t,y,y0(11),growthOn,diffusionOn),...
  'RelTol',1e-3,'AbsTol',abstol,'InitialStep',1e-3);
[t,y]= ode15s(@(t,y)singlecell.model_weberWrapper(t,y,y0(11),growthOn,diffusionOn),t,y0,opt);

