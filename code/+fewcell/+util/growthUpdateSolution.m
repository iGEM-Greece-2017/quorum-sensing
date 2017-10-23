function params= growthUpdateSolution(i,params,u)
% Update geometry parameters and init condition
% Geometry:
% - Grow new bacteria: increase nRings
% - Remove old bacteria: move bactCenter0
% Init condition:
% - Dilute old, birth new bacteria (result{2} -> params.solve.y0)
  
  % Geometry
  nBactPrev= params.g.nRings * params.g.nLayers;
  rPrevSel= params.g.startRingIdx:params.g.startRingIdx+params.g.nRings-1;
  global growth;
  params.g.nRings= params.g.nRings + growth.dr(i-1);
  ringstoDie= max(params.g.nRings - growth.maxRings, 0);
  
  % Init conditions
  nBact= params.g.nRings * params.g.nLayers;
  rNewSel= params.g.startRingIdx:params.g.startRingIdx+params.g.nRings-1;
  % Calculate dilution coefficient
  n1= sum(growth.bactRingDensity(rPrevSel)*growth.nLayers);
  n2= sum(growth.bactRingDensity(rNewSel)*growth.nLayers);
  density= repmat(growth.bactRingDensity(rPrevSel)',growth.nLayers,1);
  density= repmat(density(:)',8,1);
  yTot= sum(density.*u,2);
  yEnd= u(:,end);
  dilutionCoeff= ((n2-n1)*yEnd(2:7)./yTot(2:7)+1);
  
  % Dilute old
  u(2:7,:)= u(2:7,:)./dilutionCoeff;
  % Grow new
  uNew(2:7,:)= repmat( yEnd(2:7)./dilutionCoeff ,1,nBact-nBactPrev);
  uNew(1,:)= growth.initDNA; uNew(8,:)= 0;
  u= [u,uNew];

  density= repmat(growth.bactRingDensity(rNewSel)',growth.nLayers,1);
  density= repmat(density(:)',8,1);
  yTotNew= sum(density.*u,2);
  if ~all(abs(yTot(2:7)-yTotNew(2:7))./yEnd(2:7) < 1e-6)
    warnmsg= sprintf('[growthUpdate]: High dilution relerr: %.2g', max(abs(yTot(2:7)-yTotNew(2:7))./yEnd(2:7)));
    warning(warnmsg);
  end
  
  % Kill old
  nBactDie= ringstoDie*growth.nLayers;
  u= u(:,nBactDie+1:end);
  if ringstoDie > 0
    % move bactCenter0 <toDie> rings ahead
    params.g.bactCenter0(1)= params.g.bactCenter0(1) + ringstoDie*params.g.ringDist*params.g.bactSize(1);
    params.g.nRings= growth.maxRings;
    params.g.startRingIdx= params.g.startRingIdx + ringstoDie;
  end
  % Copy to params
  params.solve.y0= u(:);
  
  global bactProdMultiplier;
  ringX= params.g.bactCenter0(1)+(0:params.g.nRings-1)'*params.g.ringDist*params.g.bactSize(1);
  [~,bactProdMultiplier]= fewcell.util.bactRingDensity(ringX,params.g,0);
  bactProdMultiplier= repmat(bactProdMultiplier',params.g.nLayers,1); bactProdMultiplier= bactProdMultiplier(:)';
end
