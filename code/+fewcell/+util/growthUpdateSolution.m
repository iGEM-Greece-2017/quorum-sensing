function u= growthUpdateSolution(u,nBact,bactNodes)
% Each time the integration is interrupted, this function updates the system state. It dilutes the existing
% bacteria, kills inner bacteria if necessary and grows new bacteria outwards
  nBactMax= size(bactNodes,2);
  yIdx1= length(u)-nBactMax*8+1;  % first y of first bacterium
  
  %% Select bacteria that will be alive in this growth step
  global growth;
  rNewInit= growth.r0+(growth.i-2)*growth.dr+1; rPrevInit= rNewInit-growth.dr;
  toDieRing= 0;
  if rNewInit > growth.maxRings
    toDieRing= rNewInit - growth.maxRings;  % these bacteria die
  end
  growth.selBactLim= [1,nBact] + toDieRing*growth.nLayers;
  
  %% Add new rings
  startI= yIdx1+(rNewInit-1)*growth.nLayers*8;        % point to first new bacterium
  prevI= startI - growth.dr*growth.nLayers*8;     % point to first dividing bacterium (<dr> rings back)
  
  for addBactI= 1:growth.dr*growth.nLayers    % for every bacterium to add
    % Visit the bacterium it spawns from and calculate the new concentrations
    rNew= rNewInit+floor((addBactI-1)/growth.nLayers); rPrev= rPrevInit+floor((addBactI-1)/growth.nLayers);
    dilutionCoeff= growth.bactRingDensity(rPrev)/sum(growth.bactRingDensity([rPrev,rNew]));
    % Old bact: keep DNA and DNA.LuxRAHL2, dilut everything else
    yOld= [u(prevI); u(prevI+1:prevI+6)*dilutionCoeff; u(prevI+7)];
    % New bact: add init DNA, no DNA.LuxRAHL2, dilute everything else
    yNew= [growth.initDNA; u(prevI+1:prevI+6)*dilutionCoeff; 0];
    u(prevI:prevI+7)= yOld;
    u(startI:startI+7)= yNew;
    % Move start/prev indices
    startI= startI+8; prevI= prevI+8;
  end
  
end
