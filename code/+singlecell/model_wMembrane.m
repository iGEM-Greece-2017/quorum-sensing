function dydt = model_wMembrane(~,y,N)

%SI Norm Factor 
SINorms=1/60;
SINormsM=1/(60e-9);
%The kinetic constants are in SI Units (1/s, 1/sM)
aRkR=0.01*SINorms; aIkI=0.025*SINorms;
PR=6.94*SINorms; PI=6.94*SINorms;
dmR=0.347*SINorms; dmI=0.347*SINorms;
kAHL=0.04*SINorms;
k1=0.01*SINormsM; k2=1*SINorms; k3=0.05*SINormsM; k4=1*SINorms; 
k5=0.05*SINormsM; k6=10*SINorms; kR=10*SINorms; kI=2.5*SINorms;
dI=0.01*SINorms; dR=0.002*SINorms; dS=0.002*SINorms; dSS=0.002*SINorms; 
dAHL=0.01*SINorms;
Mperm=100; kr=k6/k5; k2DNA=0.015;
bactVol= 2*1*1*1e-18; spaceVol= 0.03*0.03*1e-3;
relVolume= bactVol/spaceVol;

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
