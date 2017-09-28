function dy= model_weber(t,y, p,growth)

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
dy(10,:)= p.dilution*p.colonyD*(y(6,:)-y(10,:)) -p.dA*y(10,:);
dy(1,:)= -p.klux_/p.Kdlux*y(8,:).*y(1,:) +p.klux_*y(9,:) +growth.divRate*(y(1,:)+y(9,:)) -growth.divRate*y(1,:);
dy(2,:)= p.alphaR*p.kR*y(1,:) +p.kR*y(9,:) -(growth.divRate+p.dmR)*y(2,:);
dy(3,:)= p.alphaI*p.kI*y(1,:) +p.kI*y(9,:) -(growth.divRate+p.dmI)*y(3,:);
dy(4,:)= p.pI*y(3,:) -(growth.divRate+p.dI)*y(4,:);
dy(5,:)= -p.k1_/p.Kd1*y(6,:).*y(5,:) +p.k1_*y(7,:)  +p.pR*y(2,:) -(growth.divRate+p.dR)*y(5,:);
dy(6,:)= p.k1_*y(7,:) -p.k1_/p.Kd1*y(6,:).*y(5,:) +p.kA*y(4,:) -p.D*(y(6,:)-y(10,:)) -(growth.divRate+p.dA)*y(6,:);
dy(7,:)= -p.k1_*y(7,:) +p.k1_/p.Kd1*y(6,:).*y(5,:) -2*p.k2_/p.Kd2*y(7,:).^2 +2*p.k2_*y(8,:) -(growth.divRate+p.dC)*y(7,:);
dy(8,:)= p.k2_/p.Kd2*y(7,:).^2 -p.k2_*y(8,:) -p.klux_/p.Kdlux*y(8,:).*y(1,:) +p.klux_*y(9,:) -(growth.divRate+p.dC2)*y(8,:);
dy(9,:)= p.klux_/p.Kdlux*y(8,:).*y(1,:) -p.klux_*y(9,:) -growth.divRate*y(9,:);
dy(11,:)= growth.dN;


% Proteins match up to mRNA numbers! Ref:
% - http://book.bionumbers.org/how-many-proteins-are-made-per-mrna-molecule/
end
