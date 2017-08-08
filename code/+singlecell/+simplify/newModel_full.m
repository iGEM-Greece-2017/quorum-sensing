function [dx,y]= newModel_full(t,x,u, k1,k2,a,order_f,order_b, fileArg)

dSS=0.002; k5=0.05; k6=10; k1_old= 0.01;
N= fileArg{1};
dx= singlecell.model_wMembrane(t,[x(1:6);0;x(7:9)],N);
% add back the influence of x8 on its producers
dx([5,6])= dx([5,6]) +k1_old*x(5)*x(6) -2*k1*(x(5)*x(6)).^order_f +k2*x(7).^order_b;
% calculate the new state variable
dx(8)= k1*(x(5)*x(6)).^order_f -k2*x(7).^order_b -(dSS+a)*x(7) -k5*x(7)*x(1) +k6*x(8);

dx= [dx(1:6); dx(8:10)];
y= [x(6);x(7)];
