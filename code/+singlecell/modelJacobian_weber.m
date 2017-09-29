function j= modelJacobian_weber(t,y, p,growth)

n= size(y,2);
j= zeros(11,11,n);
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
%                    *** Assumes constant p! ***
j(10,[6,10],:)= [p.dilution*p.colonyD, -p.dilution*p.colonyD*-p.dA];
j(1,[1,8,9],:)= [-p.klux_/p.Kdlux*y(8,:), -p.klux_/p.Kdlux*y(1,:), p.klux_+growth.divRate];
j(2,[1,2,9],:)= [p.alphaR*p.kR, -(growth.divRate+p.dmR), p.kR];
j(3,[1,3,9],:)= [p.alphaI*p.kI, -(growth.divRate+p.dmI), p.kI];
j(4,[3,4],:)= [p.pI, -(growth.divRate+p.dI)];
j(5,[2,5,6,7],:)= [p.pR, -p.k1_/p.Kd1*y(6,:)-(growth.divRate+p.dR), -p.k1_/p.Kd1*y(5,:), p.k1_];
j(6,[4,5,6,7,10],:)= [p.kA, -p.k1_/p.Kd1*y(6,:), -p.k1_/p.Kd1*y(5,:)-p.D-(growth.divRate+p.dA), p.k1_, p.D];
j(7,[5,6,7,8],:)= [p.k1_/p.Kd1*y(6,:), p.k1_/p.Kd1*y(5,:), -p.k1_-4*p.k2_/p.Kd2*y(7,:)-(growth.divRate+p.dC), 2*p.k2_];
j(8,[1,7,8,9],:)= [-p.klux_/p.Kdlux*y(8,:), 2*p.k2_/p.Kd2*y(7,:), -p.k2_-p.klux_/p.Kdlux*y(1,:)-(growth.divRate+p.dC2), p.klux_];
j(9,[1,8,9],:)= [p.klux_/p.Kdlux*y(8,:), p.klux_/p.Kdlux*y(1,:), -p.klux_-growth.divRate];
end
