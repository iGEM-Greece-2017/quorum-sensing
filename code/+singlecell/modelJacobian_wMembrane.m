function j= modelJacobian_wMembrane(~,y,N0, growthOn)

[aRkR,aIkI,PR,PI,dmR,dmI,kAHL, ...
 k1,k2,k3,k4,k5,k6,kR,kI,      ...
 dI,dR,dS,dSS,dAHL,Mperm,      ...
 cellDNA,growth,dilution,      ...
 diffusiveLoss,totMembPerm]= singlecell.modelCoeffs(y,N0);

% Growth model calcs
growth.f= 0;
if growthOn
  growth.f= -growth.m/growth.Nmax^growth.m * y(11).^(growth.m-1) * growth.p2;
  growth.f= growth.f + growth.p1 * growth.n * (growth.Nmin^growth.n ./ y(11).^(growth.n+1));
  growth.j= growth.r*growth.p1*growth.p2 + growth.r*y(11).*growth.f;
  growth.j= growth.j/60;    % time multiplier (hour->min)
end

j= zeros(11);
j(1, [1 8 9 11])= [-k5*y(8), -k5*y(1), k6, cellDNA*growth.j];
j(2, [1 2 9])= [aRkR, -dmR, kR];
j(3, [1 3 9])= [aIkI, -dmI, kI];
j(4, [3 4])=   [PI, -dI];
j(5, [2 5 6 7])= [PR, -k1*y(6)-dR, -k1*y(5), k2];
j(6, [4 5 6 7 10 11])= [kAHL, -k1*y(6), -k1*y(5)-dAHL-totMembPerm, k2, totMembPerm, -Mperm*(y(6)-y(10))];
j(7, [5 6 7 8])= [k1*y(6), k1*y(5), -k2-dS-4*k3*y(7), 2*k4];
j(8, [1 7 8 9])= [-k5*y(8), 2*k3*y(7), -k4-dSS-k5*y(1), k6];
j(9, [1 8 9])= [k5*y(8), k5*y(1), -k6];
j(10,[6 10 11])=  [totMembPerm*dilution, -totMembPerm*dilution-dAHL-diffusiveLoss, Mperm*dilution*(y(6)-y(10))];
j(11,11)= growth.j;

%{
global modelJacobian_minsvd;
minsvd= min(svd(j));
stability= minsvd./max(svd(j));
modelJacobian_minsvd.svd= [modelJacobian_minsvd.svd; minsvd];
modelJacobian_minsvd.stability= [modelJacobian_minsvd.stability; stability];
modelJacobian_minsvd.j= j;
%}
