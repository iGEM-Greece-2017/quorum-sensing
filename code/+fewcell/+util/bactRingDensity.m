function [n,rCoeff]= bactRingDensity(r,bactSize,lateralSpacing)
% How many bacteria of width <bactSize(1)> + <lateralSpacing> exist in a ring at <r> 
% (due to cylindrical symmetry)

% dr= bactSize(1)/2;
% a_ring= pi*((r+dr)^2 - (r-dr)^2)= 2pi*r*bactSize(1)
% a_bact= pi*dr^2 + latSpace/(2pi*r)*a_ring= pi*dr^2 + latSpace*bactSize(1)
% n= a_ring/a_bact
rCoeff= 2*pi ./ (pi/4*bactSize(1) + lateralSpacing);
n= r*rCoeff;
n= max(n, ones(length(r),1));
