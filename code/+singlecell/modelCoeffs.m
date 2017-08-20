function [aRkR,aIkI,PR,PI,dmR,dmI,kAHL, ...
          k1,k2,k3,k4,k5,k6,kR,kI,      ...
          dI,dR,dS,dSS,dAHL,Mperm,      ...
          cellDNA,growth,dilution,      ...
          diffusiveLoss,totMembPerm]= modelCoeffs(y,N0)

% The kinetic constants are in (1/min, 1/(min*nM))
% Space is in mm
aRkR=0.01; aIkI=0.025;
PR=6.94; PI=6.94;
dmR=0.347; dmI=0.347;
kAHL=0.04;
k1=0.01; k2=1; k3=0.05; k4=1; 
k5=0.05; k6=10; kR=10; kI=2.5;
dI=0.01; dR=0.002; dS=0.002; dSS=0.002; 
dAHL=0.01;
Mperm=100; kr=k6/k5; k2DNA=0.015;
cellDNA= 1.6;   % amount of DNA per cell
% Growth model constants (logistic++)
growth.r= 2.40; growth.m= 0.52; growth.n= 3.5;
growth.Nmax= 10^(9); growth.Nmin= N0*(1-1e-6);

bactVol= 2*1*1*1e-9;
% how closely together the bacteria are packed (1.05 -> separation: 5% of length per bacterium)
bactPackingFactor= 1.05;
spaceVol= bactVol*bactPackingFactor^3*y(11);
dilution= bactVol/spaceVol;
diffusiveLoss= 1e-3; % [1/min]

growth.p1= 1-(y(11)./growth.Nmax).^growth.m;
growth.p2= 1-(growth.Nmin./y(11)).^growth.n;
totMembPerm= y(11)*Mperm;
