function [x,y]= ellipseFun(x0,r,args)
% Create a horizontal ellipse around x0 with radii r. 2 line segments.

switch length(args)
  case 0
    x= 2;           % Line segments
  case 1
    bs= args{1};
    A= [0, pi;      % Parameter start vals
        pi, 2*pi;   % Parameter end vals
        1, 1;       % Left region label
        0, 0];      % Right region label
    x= A(:,bs);
  case 2
    t= args{2};
    x= x0(1)+ r(1)*cos(t);
    y= x0(2)+ r(2)*sin(t);
end
