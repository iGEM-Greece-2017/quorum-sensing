function n= bactRingDensity(r,bactSize)
% How many bacteria of width <bactSize(1)> exist in a ring at <r> 
% (due to cylindrical symmetry)

dr= bactSize(1)/2;
% pi*((r+dr)^2 - (r-dr)^2) =>
n= max(pi*4*r*dr./bactSize(1).^2, ones(length(r),1));
