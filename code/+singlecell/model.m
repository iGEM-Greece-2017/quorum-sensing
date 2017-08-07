function dydt = model(~,y)

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
        k1*y(5).*y(6) -k2*y(7) -2*k3*y(7).^2 +k4*y(8) -dS*y(7);   %7
        k3*y(7).^2 -k4*y(8) -dSS*y(8) -k5*y(8).*y(1) +k6*y(9);  %8
        k5*y(8).*y(1)-k6*y(9)];                               %9
end

%{
% Notes: LuxRAHL simplification
(y5,y6) (f1)-> y8
(y5,y6) (f2)-> y8

f1:
- y5,y6 <-> y7
- y7,y7 <-> y8
dydt= [k1*y(5).*y(6)-k2*y(7)-2*k3*y(7).^2+k4*y(8)-dS*y(7);
       2*k3*y(7).^2-k4*y(8)-dSS*y(8)];
f2:
- y5,y6 <-> y8
dydt= knew1*y(5).*y(6) - knew2*y(8) - dSS*y(8)
%}

y8'= 2k3 y7^2 -k4y8 -dss y8
y7'= k1y5y6 -k2y7 -2k3y7^2 +k4y8 -ds y7 =>
sy7-y7_0= k1y5**y6 -k2y7 -2k3y7**y7 +k4y8 -ds y7

sy8-y8_0 + sy7-y7_0= k1y5**y6 -k2y7 +k4y8 -ds y7 -k4y8 -dss y8 =>
(s+k2+ds)y7= k1y5**y6 -dssy8 -sy8+y8_0 +y7_0 =>
y7= (A) / (B) =>
sy8-y8_0= 2k3y7**y7 -k4y8 -dssy8 =>
sy8-y8_0= 2k3 ((A/B)**(A/B)) -k4y8 -dssy8 =>
y8'= 2k3 (A/B)^2 -k4y8 -dssy8

A/B= (k1y5**y6 -dssy8 -sy8+y8_0 +y7_0) / (s+k2+ds)  =>{L-1}=> 
A/B= (k1y5y6 -dssy8 -y8' +y7) ** (e^(-(k2+ds)t)*u(t))

%
y8'= 2k3 y7^2 -k4y8 -dss y8
y7'= k1y5y6 -k2y7 -2k3y7^2 +k4y8 -ds y7

y8'+y7'= +k1y5y6 -(k2+ds)y7 -dssy8






