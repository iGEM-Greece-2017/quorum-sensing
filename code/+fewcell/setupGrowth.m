function params= setupGrowth(bactRingDensity,tlist,params)
  global growth;
  global enableSinglecellEq;
  global enableGraphics;
  if params.growth.on, fprintf('\tGrowth ON\n'); else, fprintf('\tGrowth OFF\n'); end
  % Calculate growth step sizes
  if ~enableSinglecellEq, params.growth.on= false; end
  growth.on= params.growth.on;
  if growth.on
    % Copy parameters to global
    growth.r0= params.growth.r0; growth.dr= params.growth.dr;
    growth.maxRings= params.growth.maxRings; growth.nLayers= params.g.nLayers;
    growth.initDNA= params.solve.y0(1); growth.bactRingDensity= bactRingDensity;

    % Sum the number of bacteria every <dr> rings (for all layers) to calculate the growth step size
    assert(~mod(params.g.max_nRings-growth.r0, growth.dr), '<nRings>-<r0> is not divisible by <dr>; increase <nRings>');
    nSteps= (params.g.max_nRings-growth.r0) ./ params.growth.dr;
    %growth.stepSize= sum(reshape(bactRingDensity(growth.r0+1:end),params.growth.dr,[]),1);
    ringsPerGrowthstep= cell(nSteps+1,1);
    ringsPerGrowthstep{1}= [1,growth.r0];
    for i= 2:size(ringsPerGrowthstep,1)
      ringsPerGrowthstep{i}= [ringsPerGrowthstep{i-1}(1) ringsPerGrowthstep{i-1}(2)+params.growth.dr];
      if diff(ringsPerGrowthstep{i})+1 > params.growth.maxRings
        ringsPerGrowthstep{i}(1)= ringsPerGrowthstep{i}(2) - params.growth.maxRings+1;
      end
    end
    growthLevels= cellfun(@(x)sum(bactRingDensity(x(1):x(2))), ringsPerGrowthstep);
    
    %{
    maxStep= floor(params.growth.maxRings/params.growth.dr);
    if params.growth.maxRings < params.g.max_nRings
      growth.stepSize= growth.stepSize - [zeros(1,maxStep),growth.stepSize(1:end-maxStep)];
    end
    %}
    % Growth curve params
    inoculumSize= growthLevels(1);
    params.growth.gc.Nmax= growthLevels(end);
    params.growth.gc.Nmin= 0;   % disregard lag phase
    % Calculate growth curve
    growth.bactN= singlecell.growthCurve(inoculumSize,tlist,params.growth.gc);
    % Quantize the growth curve and calculate the time at which each growth step should occur
    [qgc,growth.tstep,adapt_dr]= fewcell.util.quantizeGC(growth.bactN,growthLevels,params.growth.min_dt);
    growth.dr= adapt_dr*growth.dr;

    % Plot growth curve vs quantized growth curve
    if params.viz.showGrowthCurve && enableGraphics
      figure(1); clf;
      plot(tlist/60,[growth.bactN,qgc]);
      title('Quantized growth curve');
      grid minor; axis tight; xlabel('t [h]'); ylabel('bacteria');
    end
  else
    growth.tstep= length(tlist);
  end
  minTstep= min(tlist(growth.tstep) - tlist([1; growth.tstep(1:end-1)]));
  if growth.on
    fprintf('Num of growth timesteps: %d\tMin timestep: %.0f\tQuantization mre: %.2f%%\n', ...
      length(growth.tstep),minTstep,mean((growth.bactN-qgc)./growth.bactN)*100);
    fprintf('Inoculum size: %.2g\n', inoculumSize);
  end
end
