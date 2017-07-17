function y= nan0(x)
% Turn NaN to 0 in input matrix. Useful for integrations

y= x;
y(isnan(y))= 0;
