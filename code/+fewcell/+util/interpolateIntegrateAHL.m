function [u,x,y,t,total_u]= interpolateIntegrateAHL(model,result,domainLim,zoomin,resolution,timeSubsampling,start_t)  
% Precompute interpolation for all times and store it
if nargin<7, start_t= 1; end

hexDomain= length(domainLim)==1;  % if domainLim is 2 nums, then the domain is rectangular
% X,Y grid
x= linspace(-domainLim(1)/zoomin,domainLim(1)/zoomin, resolution);
if hexDomain, y= x*sin(pi/3);
else, y= linspace(0,-domainLim(2)/zoomin,resolution);
end
[X,Y]= meshgrid(x,y);

% Slice TimeDependentResults and interpolate + integrate each slice
t= start_t: timeSubsampling: length(result.SolutionTimes);
u= zeros([size(X),length(t)]);
total_u= zeros([1,length(t)]);
i= 1;
for time= t
  resultSlice= util.sliceResult(model,result,time);
  ufun= @(x,y) util.nan0(reshape(interpolateSolution(resultSlice, x,y), size(x)));
  u(:,:,i)= ufun(X,Y);
  if hexDomain  % hexagonal domain (x,y) & CART coords
    total_u(i)= integral2(@(x,y) ufun(x,y),...
      -domainLim,domainLim, @(x)util.hexagonPerim(x,domainLim,-1), @(x)util.hexagonPerim(x,domainLim,1),...
      'AbsTol',1e-18,'RelTol',1e-6);
  else  % rectangular domain (r,z) & CYLINDRICAL coords
    total_u(i)= integral2(@(x,y) pi.*abs(x).*ufun(x,y),...
      -domainLim(1),domainLim(1), -domainLim(2),0, ...
      'AbsTol',1e-10, 'RelTol',1e-5);
  end
  i= i+1;
end
t= result.SolutionTimes(t);
