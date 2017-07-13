function [u,x,y,t,total_u]= interpolateIntegrateAHL(model,result,domainLim,zoomin,resolution,timeSubsampling,start_t)  
% Precompute interpolation for all times and store it
if nargin<7, start_t= 1; end

% X,Y grid
x= linspace(-domainLim/zoomin,domainLim/zoomin, resolution);
y= x*sin(pi/3);
[X,Y]= meshgrid(x,y);

% Slice TimeDependentResults and interpolate + integrate each slice
t= start_t: timeSubsampling: length(result.SolutionTimes);
u= zeros([size(X),length(t)]);
total_u= zeros([1,length(t)]);
i= 1;
for time= t
  resultSlice= util.sliceResult(model,result,time);
  ufun= @(x,y) reshape(interpolateSolution(resultSlice, x,y), size(x));
  u(:,:,i)= ufun(X,Y);
  total_u(i)= integral2(@(x,y) ufun(x,y),...
    -domainLim,domainLim, @(x)util.hexagonPerim(x,domainLim,-1), @(x)util.hexagonPerim(x,domainLim,1),...
    'AbsTol',1e-22,'RelTol',1e-6);
  i= i+1;
end
t= result.SolutionTimes(t);
