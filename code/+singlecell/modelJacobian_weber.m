function j= modelJacobian_weber(t,y, p,growth)

n= size(y,2);
rep_= ones(1,1,n);
y= reshape(y,11,1,n);
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
j(10,[6,10],:)= [p.dilution*p.colonyD*rep_, -p.dilution*p.colonyD*-p.dA*rep_];
j(1,[1,8,9],:)= [-p.klux_/p.Kdlux*y(8,1,:), -p.klux_/p.Kdlux*y(1,1,:), p.klux_+growth.divRate*rep_];
j(2,[1,2,9],:)= [p.alphaR*p.kR*rep_, -(growth.divRate+p.dmR)*rep_, p.kR*rep_];
j(3,[1,3,9],:)= [p.alphaI*p.kI*rep_, -(growth.divRate+p.dmI)*rep_, p.kI*rep_];
j(4,[3,4],:)= [p.pI*rep_, -(growth.divRate+p.dI)*rep_];
j(5,[2,5,6,7],:)= [p.pR*rep_, -p.k1_/p.Kd1*y(6,1,:)-(growth.divRate+p.dR), -p.k1_/p.Kd1*y(5,1,:), p.k1_*rep_];
j(6,[4,5,6,7,10],:)= [p.kA*rep_, -p.k1_/p.Kd1*y(6,1,:), -p.k1_/p.Kd1*y(5,1,:)-p.D-(growth.divRate+p.dA), p.k1_*rep_, p.D*rep_];
j(7,[5,6,7,8],:)= [p.k1_/p.Kd1*y(6,1,:), p.k1_/p.Kd1*y(5,1,:), -p.k1_-4*p.k2_/p.Kd2*y(7,1,:)-(growth.divRate+p.dC), 2*p.k2_*rep_];
j(8,[1,7,8,9],:)= [-p.klux_/p.Kdlux*y(8,1,:), 2*p.k2_/p.Kd2*y(7,1,:), -p.k2_-p.klux_/p.Kdlux*y(1,1,:)-(growth.divRate+p.dC2), p.klux_*rep_];
j(9,[1,8,9],:)= [p.klux_/p.Kdlux*y(8,1,:), p.klux_/p.Kdlux*y(1,1,:), -p.klux_-growth.divRate*rep_];
end
