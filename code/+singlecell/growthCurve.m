function n= growthCurve(n0,bactNmax,tlist)

  % Growth parameters
  growth.r= 2.2/60;
  growth.m= 0.52; growth.n= 3.5;
  growth.Nmax= bactNmax;
  %growth.Nmin= n0*(1-1e-6);
  growth.Nmin= 0;   % disregard lag phase
  
  opt= odeset('RelTol',1e-5,'InitialStep',1e-3);
  [~,n]= ode15s(@(t,y)growthWrapper(t,y,growth),tlist,n0,opt);
end

function dN= growthWrapper(~,n,growth)
  % Growth calculations
  p1= 1-(n./growth.Nmax).^growth.m;
  p2= 1-(growth.Nmin./n).^growth.n;
  dN= growth.r*n*p1*p2;
end
