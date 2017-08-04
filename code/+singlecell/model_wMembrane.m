function dydt = model_wMembrane(~,y,N)

% The kinetic constants are in (1/min, 1/(min*nM))
% Space is in μm
aRkR=0.01; aIkI=0.025;
PR=6.94; PI=6.94;
dmR=0.347; dmI=0.347;
kAHL=0.04;
k1=0.01; k2=1; k3=0.05; k4=1; 
k5=0.05; k6=10; kR=10; kI=2.5;
dI=0.01; dR=0.002; dS=0.002; dSS=0.002; 
dAHL=0.01;
Mperm=100; kr=k6/k5; k2DNA=0.015;
%bactVol= 2*1*1; spaceVol= 10*5*5;
%relVolume= bactVol/spaceVol;

relVolume= 1;

%y(1)=DNA
%y(2)=mRNALuxR
%y(3)=mRNALuxI
%y(4)=LuxI
%y(5)=LuxR
%y(6)=AHLint
%y(7)=LuxRAHL
%y(8)=LuxRAHL2
%y(9)=DNALuxRAHL2
%y(10)=AHLext
dydt = [-k5*y(8)*y(1)+k6*y(9);
        aRkR*y(1)-dmR*y(2)+kR*y(9);
        aIkI*y(1)-dmI*y(3)+kI*y(9);
        PI*y(3)-dI*y(4);
        PR*y(2)-k1*y(5)*y(6)+k2*y(7)-dR*y(5); %5
        kAHL*y(4)-k1*y(5)*y(6)+k2*y(7)-dAHL*y(6)-Mperm*(y(6)-y(10));
        k1*y(5)*y(6)-k2*y(7)-2*k3*y(7)^2+k4*y(8)-dS*y(7);
        2*k3*y(7)^2-k4*y(8)-dSS*y(8)-k5*y(8)*y(1)+k6*y(9);
        k5*y(8)*y(1)-k6*y(9);
        relVolume*N*Mperm*(y(6)-y(10))-dAHL*y(10)];
end
