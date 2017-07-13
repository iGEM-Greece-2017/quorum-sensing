function [x,y]= hexFun(r,args)
% Create a horizontal regular hexagon around (0,0) with radius r. 6 line segments.

A= [-2, -1.5, -0.5, 0, 0.5, 1.5;   % Parameter start vals
    -1.5, -0.5, 0, 0.5, 1.5, 2;    % Parameter end vals
    1, 1, 1, 1, 1, 1;       % Left region label
    0, 0, 0, 0, 0, 0];      % Right region label
switch length(args)
  case 0
    x= 6;           % Line segments
  case 1
    bs= args{1};
    x= A(:,bs);
  case 2
    %bs= args{1}
    t= args{2};
    %size(t)
    x= zeros(size(t));
    y= zeros(size(t));
    x(t<=0)= r*(-t(t<=0)-1);
    x(t>0)= r*(t(t>0)-1);
    
    % Define 6 areas based on t
    area{1}= t>=A(1,1) & t<A(2,1);
    area{2}= t>=A(1,2) & t<A(2,2);
    area{3}= t>=A(1,3) & t<A(2,3);
    area{4}= t>=A(1,4) & t<A(2,4);
    area{5}= t>=A(1,5) & t<A(2,5);
    area{6}= t>=A(1,6) & t<=A(2,6);
    % Assign y
    y(area{1})= r*sin(pi/3)/(r*cos(pi/3)-r) * (x(area{1})-r);
    y(area{2})= r*sin(pi/3);
    y(area{3})= r*sin(pi/3)/(r-r*cos(pi/3)) * (x(area{3})+r);
    y(area{4})= -(r*sin(pi/3)/(r-r*cos(pi/3)) * (x(area{4})+r));
    y(area{5})= -r*sin(pi/3);
    y(area{6})= -(r*sin(pi/3)/(r*cos(pi/3)-r) * (x(area{6})-r));
end
