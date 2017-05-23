%% Parameters
params= struct(...
  'dt',4e-3,...               % Base world timestep. Corresponds to diffusion timescale
  'r',1e-6,...                % Hex grid cell radius
  'initDNA',1.6E-7...        % Initial [DNA] in each bacterium
  );
params.period= struct(...     % How many times slower these processes are than base timestep
  'singleBactModel',250,...
  'replication',250*60*20 ...
  );
params.diffusion= struct(...
  'D_AHL',4.9E-10, ...        % [m2/s]
  'C_agar',0.9 ...            % [no unit]
  );
worldSize= [60 60];
world= multicell.World(worldSize,0.01,params);

time_stop= 0.05;   % Time to stop [min]
%% Simulate world
steps= time_stop*60/params.dt;
hist= struct('hasBact',cell(steps,1),'AHL',0,'bactState',0);
for i=1:steps
  [~,~,hist(i).hasBact,hist(i).AHL,hist(i).bactState]= world.getState();
  world.step();
end

%% Visualize results
% Total AHL plot
totalAHL= arrayfun(@(x) sum(x.AHL(:)),hist);
figure(1);
plot(totalAHL);
pause(1);

% AHL concentration heatmap movie
figure(2);
for i=1:steps
  ahl= hist(i).AHL;
  ahl_g= imagesc(ahl, [min(ahl(:)),max(ahl(:))+10*eps]);
  title(num2str(i));
  drawnow;
  %pause(0.005);
end
