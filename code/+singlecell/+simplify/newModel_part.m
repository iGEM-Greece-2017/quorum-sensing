function [dx,y]= newModel_part(t,x,u, k1,k2,a,order,~)

dSS=0.002; k5=0.05; k6=10;

% u(1,2) <- y(5)*y(6)
% u(3) <- y(1)
% u(4) <- y(9)
% x <- y(8)

dx= k1*(u(1).*u(2)).^order -k2*x -(dSS+a)*x -k5*x*u(3) +k6*u(4);
y= [x; (-2*k1*(u(1).*u(2)).^order +k2*x)];   % #2: amount of reagents (y5,y6) consumed
