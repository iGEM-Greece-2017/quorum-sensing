function total_u= integrate(ufun,domainLim,absTol)

total_u= integral2(@(x,y) pi.*abs(x).*ufun(x,y),...
  -domainLim(1),domainLim(1), -domainLim(2),0, ...
  'AbsTol',absTol, 'RelTol',1e-3);