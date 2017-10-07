% Display precalculated result
%load('data/tmpresults.mat');

%% Plot solution
global enableGraphics; enableGraphics= true;
% Prepare solution interpolation
fprintf('--> Interpolating solution...\n');
[AHLDistrib,x,y,interp_t,totalAHL]= ...
  fewcell.util.interpolateIntegrateAHL(model,result,params,470,714);
% Plot
if params.viz.figID(1), fprintf('--> [paused] Press key to plot the solution\n'); pause; end
fprintf('--> Plotting solution...\n');
fewcell.output(model,result,tlist,AHLDistrib,x,y,interp_t,totalAHL,params.viz);
