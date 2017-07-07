function y= hexagonPerim(x,r,oneSide)
% Returns the up and down y values of a regular hexagon centered around 0 with radius r 
% and 0 angle with x-axis
%
% Input:
% - x {N}: x coordinate vector
% - r: radius
% - oneSide: if exists, output only 1 side. 1->upper, -1->lower
% Output:
% - y {N}x2: up and down value of the hexagon's perimeter

% Cut the range of x into 5 areas
area{1}= x < -r;
area{2}= x >= -r & x < -cos(pi/3)*r;
area{3}= x >= -cos(pi/3)*r & x < cos(pi/3)*r;
area{4}= x >= cos(pi/3)*r & x < r;
area{5}= x >= r;
% Assign the hexagon's perimeter values for the upper side
y= zeros(length(x),2);
y(area{1},1)= 0;
y(area{2},1)= r*sin(pi/3)/(r-r*cos(pi/3)) * (x(area{2})+r);
y(area{3},1)= r*sin(pi/3);
y(area{4},1)= r*sin(pi/3)/(r*cos(pi/3)-r) * (x(area{4})-r);
y(area{5},1)= 0;
% Lower side
y(:,2)= -y(:,1);
% Check if only 1 side requested
if nargin==3
  if oneSide==1, y= y(:,1);
  elseif oneSide==-1, y= y(:,2);
  end
end
% Same dimension as x
if size(x,1)==1, y= y'; end
