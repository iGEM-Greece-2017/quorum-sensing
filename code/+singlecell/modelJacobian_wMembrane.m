function j= modelJacobian_wMembrane(~,y,N)

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


j= zeros(9);
j(1, [1 8 9])= [-k5*y(8), -k5*y(1), k6];
j(2, [1 2 9])= [aRkR, -dmR, kR];
j(3, [1 3 9])= [aIkI, -dmI, kI];
j(4, [3 4])=   [PI, -dI];
j(5, [2 5 6 7])= [PR, -k1*y(6)-dR, -k1*y(5), k2];
j(6, [4 5 6 7 10])= [kAHL, -k1*y(6), -k1*y(5)-dAHL-Mperm, k2, Mperm];
j(7, [5 6 7 8])= [k1*y(6), k1*y(5), -k2-dS-4*k3*y(7), k4];
j(8, [1 7 8 9])= [-k5*y(8), 4*k3*y(7), -k4-dSS-k5*y(1), k6];
j(9, [1 8 9])= [k5*y(8), k5*y(1), -k6];
j(10,[6 10])=  [N*Mperm*relVolume, -N*Mperm*relVolume-dAHL];

