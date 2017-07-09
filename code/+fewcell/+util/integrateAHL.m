function total_u= integrateAHL(result,domainLim,start_t)
% Precompute interpolation and integral for all times and store it
ufun= @(x,y,t) reshape(interpolateSolution(result, x,y,t), size(x));
total_u= zeros(1,length(result.SolutionTimes));
if nargin<4, start_t= 1; end

for ti= start_t:length(result.SolutionTimes)
  total_u(ti)= integral2(@(x,y) ufun(x,y,ti),...
    -domainLim,domainLim, @(x)util.hexagonPerim(x,domainLim,-1), @(x)util.hexagonPerim(x,domainLim,1),...
    'AbsTol',0,'RelTol',4);
end
