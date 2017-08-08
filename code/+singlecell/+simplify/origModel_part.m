function [dx,y]= origModel_part(t,x,u,~)

k1=0.01; k2=1; k3=0.05; k4=1; k5=0.05; k6=10;
dS=0.002; dSS=0.002; 

% u(1,2) <- y(5)*y(6)
% u(3) <- y(1)
% u(4) <- y(9)
% x(1) <- y(7)
% x(2) <- y(8)

dx= [(k1*u(1).*u(2) -k2*x(1) -2*k3*x(1).^2 +k4*x(2) -dS*x(1));    %7
     (k3*x(1).^2 -k4*x(2) -dSS*x(2) -k5*x(2)*u(3) +k6*u(4))];     %8
y= [x(2); (-k1*u(1).*u(2) +k2*x(1))];   % #2: amount of reagents (y5,y6) consumed
