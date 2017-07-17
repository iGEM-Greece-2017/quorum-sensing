function [x,y]= rectangleFun(upperLeft,l,args)
% Create a horizontal rectangle starting at the <upperLeft> corner
% and with sides of length <l(1,2)>

switch length(args)
  case 0
    x= 4;
  case 1
    bs= args{1};
    A= [0,   l(1),  sum(l),     sum(l)+l(1);  % Parameter start vals
        l(1),sum(l),sum(l)+l(1),2*sum(l);     % Parameter end vals
        0, 0, 0, 0;       % Left region label
        1, 1, 1, 1];      % Right region label
    x= A(:,bs);
  case 2
    t= args{2};
    % upper side
    area{1}= t<=l(1);
    x(area{1})= upperLeft(1) + t(area{1});
    y(area{1})= upperLeft(2);
    % right side
    area{2}= t>l(1) & t<=sum(l);
    x(area{2})= upperLeft(1) + l(1);
    y(area{2})= upperLeft(2) - (t(area{2})-l(1));
    % lower side
    area{3}= t>sum(l) & t<=sum(l)+l(1);
    x(area{3})= upperLeft(1) + l(1) - (t(area{3})-sum(l));
    y(area{3})= upperLeft(2) - l(2);
    % left side
    area{4}= t>sum(l)+l(1);
    x(area{4})= upperLeft(1);
    y(area{4})= upperLeft(2) - l(2) + (t(area{4})-sum(l)-l(1));    
end
