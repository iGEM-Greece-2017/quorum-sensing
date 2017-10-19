function [params,tlist,bactRingDensity_allRings]= initSetup(params)
  tlist= linspace(params.t.tstart, params.t.tstop, params.t.timePoints);
  % Update parameters
  params.g.nRings= params.growth.r0;
  params.g.bactCenter0= params.g.init_bactCenter0;
  params.g.startRingIdx= 1;
  % Calculate the number of simulated bacteria on each ring
  global bactProdMultiplier;
  ringX= params.g.init_bactCenter0(1)+(0:params.g.max_nRings-1)'*params.g.ringDist*params.g.bactSize(1);
  [bactRingDensity_allRings,bactProdMultiplier]= fewcell.util.bactRingDensity(ringX,params.g,0);
  % simulate a higher density, counteracting the spacing between bacteria
  bactRingDensity_allRings= bactRingDensity_allRings*params.g.nLayers.*bactProdMultiplier;
  totalBacteria= sum(bactRingDensity_allRings);
  if params.growth.on && params.growth.maxRings < params.g.nRings
    totalBacteria= sum(bactRingDensity_allRings(end-params.growth.maxRings:end));
  end
  fprintf('Bacteria: %.3g', round(totalBacteria));
  
  params= fewcell.setupGrowth(bactRingDensity_allRings,tlist,params);
  
  ringX= params.g.bactCenter0(1)+(0:params.g.nRings-1)'*params.g.ringDist*params.g.bactSize(1);
  [~,bactProdMultiplier]= fewcell.util.bactRingDensity(ringX,params.g,0);
  bactProdMultiplier= repmat(bactProdMultiplier',params.g.nLayers,1); bactProdMultiplier= bactProdMultiplier(:)';
end
