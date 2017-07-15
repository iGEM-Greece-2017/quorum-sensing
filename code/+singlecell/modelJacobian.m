function j= modelJacobian(~,y)

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

j= zeros(9);
j(1, [1 8 9])= [-k5*y(8), -k5*y(1), k6];
j(2, [1 2 9])= [aRkR, -dmR, kR];
j(3, [1 3 9])= [aIkI, -dmI, kI];
j(4, [3 4])=   [PI, -dI];
j(5, [2 5 6 7])= [PR, -k1*y(6)-dR, -k1*y(5), k2];
j(6, [4 5 6 7])= [kAHL, -k1*y(6), -k1*y(5)-dAHL, k2];
j(7, [5 6 7 8])= [k1*y(6), k1*y(5), -k2-dS-4*k3*y(7), k4];
j(8, [1 7 8 9])= [-k5*y(8), 4*k3*y(7), -k4-dSS-k5*y(1), k6];
j(9, [1 8 9])= [k5*y(8), k5*y(1), -k6];

end
