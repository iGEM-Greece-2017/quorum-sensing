function [u,nBact]= growthUpdateSolution(i,nBact,result)
% Each time the integration is interrupted, this function updates the system state. It dilutes the existing
% bacteria, kills inner bacteria if necessary and grows new bacteria outwards
  global bactNodes;
  global yyResults;
  global growth;
  if ~isempty(result)
    u= result.NodalSolution(end,:);
  else
    u= 0;
  end
  nBactMax= size(bactNodes,2);
  yIdx1= length(u)-nBactMax*8+1;  % first y of first bacterium
  rNewInit= growth.r0+(i-2)*growth.dr+1;
  
  %% Select bacteria that will be alive in this growth step
  % how many rings will die
  toDieRing= 0;
  if rNewInit > growth.maxRings
    toDieRing= rNewInit - growth.maxRings;  % these bacteria die
  end
  nBact= nBact - toDieRing*growth.nLayers;
  growth.selBactLim= [1,nBact]+toDieRing*growth.nLayers;
  
  %% Dilute old bacteria
  % Dilution protocol (to simulate cytoplasm growth):
  % Assume all bacteria participate in each growth step, not just the outermost layers
  % The total amount of each species must be maintained during this redistribution
  % Calculate a dilution factor <d> to dilute all the existing bacteria.
  % The concentrations in the new bacteria will be the same as those of the previous outermost ring
  % after dilution
  % 
  % Total bacteria before/after
  rPrevSel= floor((growth.selBactLim+1)/growth.nLayers); rPrevSel(2)= rPrevSel(2) - growth.dr;
  rPrevSel= rPrevSel(1):rPrevSel(2); % previously selected rings:
  nNewBact= growth.dr*growth.nLayers; nOldBact= nBact-nNewBact;
  n1= sum(growth.bactRingDensity(rPrevSel)*growth.nLayers);
  n2= n1+sum(growth.bactRingDensity(rNewInit:rNewInit+growth.dr-1)*growth.nLayers);
  % Total amounts of each y (analogous to sum of concentrations)
  bactIdx= yIdx1+(growth.selBactLim(1)-1)*8; bactIdx= [bactIdx, bactIdx+(nOldBact-1)*8];
  density= repmat(growth.bactRingDensity(rPrevSel),1,growth.nLayers)';
  density= repmat(density(:)',8,1);
  yTot= sum(density .* reshape( u(bactIdx(1):bactIdx(2)+7) ,8,[]),2);  % TODO: growth.selBactLim(2) -> teleftaio i meta
  yEnd= u(bactIdx(2): bactIdx(2)+7);
  dilutionCoeff= ((n2-n1)*yEnd(2:7)./yTot(2:7)+1).^(-1);
  % Only dilute y2-y8 (not DNA or DNA.LuxRAHL2)
  yDiluteIdx= reshape(bactIdx(1):bactIdx(2)+7, 8,[]);
  yDiluteIdx= yDiluteIdx(2:7,:); yDiluteIdx= yDiluteIdx(:);
  u(yDiluteIdx)= repmat(dilutionCoeff,nOldBact,1).*u(yDiluteIdx);
  
  %% Grow new bacteria
  newIdx= bactIdx(2)+8; newIdx= [newIdx, newIdx+(nNewBact-1)*8];
  yDiluteIdx= reshape(newIdx(1):newIdx(2)+7, 8,[]);
  yDiluteIdx= yDiluteIdx(2:7,:); yDiluteIdx= yDiluteIdx(:);
  u(yDiluteIdx)= repmat(dilutionCoeff.*yEnd(2:7), nNewBact,1);
  % Initialize their DNA (DNA.LuxRAHL2 is already 0)
  yDNAIdx= reshape(newIdx(1):newIdx(2)+7, 8,[]);
  yDNAIdx= yDNAIdx(1,:); yDNAIdx= yDNAIdx(:);
  u(yDNAIdx)= growth.initDNA;
  
      densNew= repmat(growth.bactRingDensity(rNewInit:rNewInit+growth.dr-1),1,growth.nLayers)';
      densNew= repmat(densNew(:)',8,1);
      yTotOld2= sum(density .* reshape( u(bactIdx(1):bactIdx(2)+7) ,8,[]),2);
      yTotNew2= sum(densNew .* reshape( u(newIdx(1):newIdx(2)+7) ,8,[]),2);
      yTot2= yTotOld2+yTotNew2;
      assert(all(abs(yTot(2:7)-yTot2(2:7))./yEnd(2:7) < 1e-8), 'Dilution error!');
  
end
