function n= growthCurve(n0,tlist, growthParams)
  opt= odeset('RelTol',1e-5,'InitialStep',1e-3);
  [~,n]= ode15s(@(t,y)growthWrapper(t,y,growthParams),tlist,n0,opt);
end

function dN= growthWrapper(~,n,growth)
  % Growth calculations
  p1= 1-(n./growth.Nmax).^growth.m;
  p2= 1-(growth.Nmin./n).^growth.n;
  dN= growth.r*n*p1*p2;
end
