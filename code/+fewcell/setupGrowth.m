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
    maxStep= floor(params.growth.maxRings/params.growth.dr);

    % Sum the number of bacteria every <dr> rings (for all layers) to calculate the growth step size
    assert(~mod(params.g.max_nRings-growth.r0, growth.dr), '<nRings>-<r0> is not divisible by <dr>; increase <nRings>');
    growth.stepSize= sum(reshape(bactRingDensity(growth.r0+1:end),params.growth.dr,[]),1);
    if params.growth.maxRings < params.g.max_nRings
      growth.stepSize= growth.stepSize - [zeros(1,maxStep),growth.stepSize(1:end-maxStep)];
    end
    % Growth curve params
    inoculumSize= sum(bactRingDensity(1:params.growth.r0));
    params.growth.gc.Nmax= sum(growth.stepSize)+inoculumSize;
    params.growth.gc.Nmin= 0;   % disregard lag phase
    % Calculate growth curve
    growth.bactN= singlecell.growthCurve(inoculumSize,tlist,params.growth.gc);
    % Calculate the time at which each growth step should occur
    [qgc,growth.tstep]= fewcell.util.quantizeGC(growth.bactN,inoculumSize,growth.stepSize);
    growth.tstep= [growth.tstep;length(tlist)];

    % Plot growth curve vs quantized growth curve
    if params.viz.showGrowthCurve && enableGraphics
      figure(1); plot(tlist/60,[growth.bactN,qgc]);
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
