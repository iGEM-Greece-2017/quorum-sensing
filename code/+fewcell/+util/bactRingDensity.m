function [n,prodMultiplier]= bactRingDensity(r,geom,lateralSpacing)
% How many bacteria of width <bactSize(1)> + <lateralSpacing> exist in a ring at <r> 
% (due to cylindrical symmetry)

rCoeff= 2*pi ./ (pi/4*geom.bactSize(1) + lateralSpacing);
n= r*rCoeff;
n= max(n, ones(length(r),1));

nNext= n(end)+n(end)-n(end-1);
dn= diff([n;nNext])/(1+geom.ringDist);
prodMultiplier= 1+geom.ringDist + geom.ringDist*(geom.ringDist+1)/2 * dn./n;
prodMultiplier= prodMultiplier*(1+geom.layerSeparation);
