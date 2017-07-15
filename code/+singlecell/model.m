function dydt = model(~,y)

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
%y(1)=DNA
%y(2)=mRNALuxR
%y(3)=mRNALuxI
%y(4)=LuxI
%y(5)=LuxR
%y(6)=AHLint
%y(7)=LuxRAHL
%y(8)=LuxRAHL2
%y(9)=DNALuxRAHL2

dydt = [-k5*y(8).*y(1)+k6*y(9);               %1
        aRkR*y(1)-dmR*y(2)+kR*y(9);           %2
        aIkI*y(1)-dmI*y(3)+kI*y(9);           %3
        PI*y(3)-dI*y(4);                      %4
        PR*y(2)-k1*y(5).*y(6)+k2*y(7)-dR*y(5);                %5
        kAHL*y(4)-k1*y(5).*y(6)+k2*y(7)-dAHL*y(6);            %6
        k1*y(5).*y(6)-k2*y(7)-2*k3*y(7).^2+k4*y(8)-dS*y(7);   %7
        2*k3*y(7).^2-k4*y(8)-dSS*y(8)-k5*y(8).*y(1)+k6*y(9);  %8
        k5*y(8).*y(1)-k6*y(9)];                               %9
end
