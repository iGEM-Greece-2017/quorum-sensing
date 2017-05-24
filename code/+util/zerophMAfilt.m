function y = zerophMAfilt(winL,x)
% Zero-phase MA filter for x with the specified window length and 0-order padding

if winL > 1
  assert(size(x,2)==1);   % assume x is a column
  padL= floor(winL/2);
  padx= [repmat(x(1),padL,1);x;repmat(x(end),padL,1)];
  idx= repmat((1:winL)',1,length(x)) + repmat(0:length(x)-1,winL,1);
  y= mean(padx(idx))';
  y= y(1:length(x));
else
  y= x;
end
end

