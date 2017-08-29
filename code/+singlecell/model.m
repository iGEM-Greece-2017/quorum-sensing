function dydt = model(~,y)

[aRkR,aIkI,PR,PI,dmR,dmI,kAHL, ...
 k1,k2,k3,k4,k5,k6,kR,kI,      ...
 dI,dR,dS,dSS,dAHL]= singlecell.modelCoeffs([y;1;1],1);

%y(1)=DNA
%y(2)=mRNALuxR
%y(3)=mRNALuxI
%y(4)=LuxI
%y(5)=LuxR
%y(6)=AHLint
%y(7)=LuxRAHL
%y(8)=LuxRAHL2
%y(9)=DNALuxRAHL2
dydt = [-k5*y(8).*y(1)+k6*y(9);                 %1
        aRkR*y(1)-dmR*y(2)+kR*y(9);             %2
        aIkI*y(1)-dmI*y(3)+kI*y(9);             %3
        PI*y(3)-dI*y(4);                        %4
        PR*y(2)-k1*y(5).*y(6)+k2*y(7)-dR*y(5);                %5
        kAHL*y(4)-k1*y(5).*y(6)+k2*y(7)-dAHL*y(6);            %6
        k1*y(5).*y(6)-k2*y(7)-2*k3*y(7).^2+2*k4*y(8)-dS*y(7); %7
        k3*y(7).^2-k4*y(8)-dSS*y(8)-k5*y(8).*y(1)+k6*y(9);    %8
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







