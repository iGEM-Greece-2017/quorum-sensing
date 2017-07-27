function x= nan0(x)
% Turn NaN & <0 to 0 in input matrix. Useful for integrations

x(isnan(x))= 0;
x(x<0)= 0;
