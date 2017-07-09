function [u,x,y]= interpolateAHL(result,domainLim,resolution,start_t)  
% Precompute interpolation and integral for all times and store it
x= linspace(-domainLim,domainLim, resolution);
y= x*sin(pi/3);
[X,Y]= meshgrid(x,y);
ufun= @(x,y,t) reshape(interpolateSolution(result, x,y,t), size(x));
u= zeros([size(X),length(result.SolutionTimes)]);
if nargin<4, start_t= 1; end

for ti= start_t:length(result.SolutionTimes)
  u(:,:,ti)= ufun(X,Y,ti);
end
