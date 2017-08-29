function dy= model_weber(t,y, growthOn, diffusionOn)

[Kd1,k1_,Kd2,k2_,kA,Kdlux,klux_,kR,kI,        ...
 pR,pI,alphaR,alphaI,dA,dC2,dC,dR,dI,dmR,dmI, ...
 D,dilution,                         ...
 growth,diffusiveLossRate,colonyD]= singlecell.modelCoeffs_weber(y,growthOn,diffusionOn);
dy= zeros(11,1);

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
%y(11)=bacteria count
dy(1)= -klux_/Kdlux*y(8)*y(1) +klux_*y(9) +growth.divRate*(y(1)+y(9)) -growth.divRate*y(1);
dy(2)= alphaR*kR*y(1) +kR*y(9) -(growth.divRate+dmR)*y(2);
dy(3)= alphaI*kI*y(1) +kI*y(9) -(growth.divRate+dmI)*y(3);
dy(4)= pI*y(3) -(growth.divRate+dI)*y(4);
dy(5)= -k1_/Kd1*y(6)*y(5) +k1_*y(7)  +pR*y(2) -(growth.divRate+dR)*y(5);
dy(6)= k1_*y(7) -k1_/Kd1*y(6)*y(5) +kA*y(4) -D*(y(6)-y(10)) -(growth.divRate+dA)*y(6);
dy(7)= -k1_*y(7) +k1_/Kd1*y(6)*y(5) -2*k2_/Kd2*y(7).^2 +2*k2_*y(8) -(growth.divRate+dC)*y(7);
dy(8)= k2_/Kd2*y(7).^2 -k2_*y(8) -klux_/Kdlux*y(8)*y(1) +klux_*y(9) -(growth.divRate+dC2)*y(8);
dy(9)= klux_/Kdlux*y(8)*y(1) -klux_*y(9) -growth.divRate*y(9);
dy(10)= dilution*colonyD*(y(6)-y(10)) -(diffusiveLossRate+dA)*y(10);
dy(11)= growth.dN;


% Proteins match up to mRNA numbers! Ref:
% - http://book.bionumbers.org/how-many-proteins-are-made-per-mrna-molecule/
end
