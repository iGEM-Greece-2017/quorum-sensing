function j= modelJacobian(~,y)

[aRkR,aIkI,PR,PI,dmR,dmI,kAHL, ...
 k1,k2,k3,k4,k5,k6,kR,kI,      ...
 dI,dR,dS,dSS,dAHL]= singlecell.modelCoeffs([y;1;1],1);

j= zeros(9);
j(1, [1 8 9])= [-k5*y(8), -k5*y(1), k6];
j(2, [1 2 9])= [aRkR, -dmR, kR];
j(3, [1 3 9])= [aIkI, -dmI, kI];
j(4, [3 4])=   [PI, -dI];
j(5, [2 5 6 7])= [PR, -k1*y(6)-dR, -k1*y(5), k2];
j(6, [4 5 6 7])= [kAHL, -k1*y(6), -k1*y(5)-dAHL, k2];
j(7, [5 6 7 8])= [k1*y(6), k1*y(5), -k2-dS-4*k3*y(7), 2*k4];
j(8, [1 7 8 9])= [-k5*y(8), 2*k3*y(7), -k4-dSS-k5*y(1), k6];
j(9, [1 8 9])= [k5*y(8), k5*y(1), -k6];

end
