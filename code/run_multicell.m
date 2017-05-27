global global_graphicsOn;
global_graphicsOn= true;
%% Parameters
params= struct(...
  'dt',7e-4,...               % Base world timestep. Corresponds to diffusion timescale
  'r',1e-6,...                % Hex grid cell radius
  'initDNA',1.6e-7 ...        % Initial [DNA] in each bacterium
  );
params.period= struct(...     % How many times slower these processes are than base timestep
  'singleBactModel',500,...
  'replication',250*60*20 ...
  );
params.diffusion= struct(...
  'D_AHL',4.9e-10, ...        % [m2/s]
  'C_agar',0.9 ...            % [no unit]
  );
worldSize= [13 13];
world= multicell.World(worldSize,2e-2,params);
%load('results/2/bactPos.mat');
world.hasBacterium= false(worldSize);
world.hasBacterium(6,6)= true; world.hasBacterium(8,8)= true;

% Visualize bacterial positions
if global_graphicsOn
  %figure(4); spy(world.hasBacterium); title('Bacterial positions');
  %pause(1);
end
tic
time_stop= 150;   % Time to stop [min]
%% Simulate world
histSample= 0.05;
steps= floor(time_stop*60/params.dt+3*eps);
%hist= struct('hasBact',cell(ceil(steps*histSample),1),'AHL',0,'bactState',0);
hist= struct('AHL',cell(ceil(steps*histSample),1));
for i=1:steps
  world.step();
  if ~mod(i,1/histSample)
    histIdx= i*histSample;
    [~,~,bactPos,AHL,~]= world.getState();
    hist(histIdx).AHL= single(AHL);
  end
end
toc
%% Visualize results
if global_graphicsOn
  % Total AHL plot
  figure(1); totAHL= viz.totalAHL(hist(1:end-1),params.period.singleBactModel*histSample);
  fprintf('[main]: max total AHL: %.4g\n', max(totAHL(:,2)));

  % Intracellular AHL plot
  %figure(2); incellAHL= viz.incellAHL(hist(1:end-1),5);

  % AHL concentration heatmap movie
  figure(3); viz.makeMovie_spatialAHL(hist(steps*histSample/2:2:end),steps*histSample/2);
end
