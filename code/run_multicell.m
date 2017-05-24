global global_graphicsOn;
global_graphicsOn= true;
%% Parameters
params= struct(...
  'dt',1e-4,...               % Base world timestep. Corresponds to diffusion timescale
  'r',1e-6,...                % Hex grid cell radius
  'initDNA',1.6E-7...        % Initial [DNA] in each bacterium
  );
params.period= struct(...     % How many times slower these processes are than base timestep
  'singleBactModel',1500,...
  'replication',250*60*20 ...
  );
params.diffusion= struct(...
  'D_AHL',4.9E-10, ...        % [m2/s]
  'C_agar',0.9 ...            % [no unit]
  );
worldSize= [100 100];
world= multicell.World(worldSize,3e-3,params);
load('results/2/bactPos.mat'); world.hasBacterium= bactPos;

% Visualize bacterial positions
if global_graphicsOn
  %figure(4); spy(world.hasBacterium); title('Bacterial positions');
  %pause(1);
end
tic
time_stop= 1;   % Time to stop [min]
%% Simulate world
steps= floor(time_stop*60/params.dt+3*eps);
hist= struct('hasBact',cell(ceil((1-0.95)*steps),1),'AHL',0,'bactState',0);
for i=1:steps
  world.step();
  if i>0.95*steps
    histIdx= ceil(i-0.95*steps);
    [~,~,hist(histIdx).hasBact,hist(histIdx).AHL,hist(histIdx).bactState]= world.getState();
  end
end
toc
%% Visualize results
if global_graphicsOn
  % Total AHL plot
  figure(1); totAHL= viz.totalAHL(hist(1:end-1),params.period.singleBactModel);
  fprintf('[main]: max total AHL: %.4g\n', max(totAHL(:,2)));

  % Intracellular AHL plot
  figure(2); incellAHL= viz.incellAHL(hist(1:end-1),5);

  % AHL concentration heatmap movie
  %figure(3); viz.makeMovie_spatialAHL(hist(steps/2:2:end),steps/2);
end
