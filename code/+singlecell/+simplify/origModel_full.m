function [dx,y]= origModel_full(t,x,u, fileArg)
N= fileArg{1};
dx= singlecell.model_wMembrane(t,x,N);
y= [x(6);x(8)];
