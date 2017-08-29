function [u,x,y,t,total_u]= interpolateIntegrateAHL(model,result, ...
            params,startTimeIdx,stopTimeIdx)
          
% Precompute interpolation for all times and store it
if nargin<4, startTimeIdx= 1; stopTimeIdx= length(result.SolutionTimes); end
if nargin==4, error('[interpolateAHL]: specify end time'); end
zoomDomainLim= params.g.domainLim ./ params.viz.zoominFactor;
resolution= params.viz.interpResolution;

% X,Y grid
if length(resolution)==1, resolution(2)= resolution(1); end
x= linspace(0,zoomDomainLim(1), resolution(1));
y= linspace(0,-zoomDomainLim(2),resolution(2));
[X,Y]= meshgrid(x,y);

% Slice TimeDependentResults and interpolate + integrate each slice
tic;
timePoints= min(params.viz.timePoints,length(result.SolutionTimes));
t= floor(linspace(startTimeIdx,stopTimeIdx, timePoints));
u= zeros([size(X),length(t)]);
total_u= zeros([1,length(t)]);
if ~params.viz.figID(1), return; end
%defaultAbsTol= 1e-1;
i= 1;
resultSlice= cell(timePoints,1);
for time= t
  resultSlice{time}= util.sliceResult(model,result,time);
end
for time= t
  ufun= @(x,y) util.nan0(reshape(interpolateSolution(resultSlice{time}, x,y), size(x)));
  u(:,:,i)= ufun(X,Y);
  % rectangular domain (r,z) & CYLINDRICAL coords
  %uSlice= u(:,:,i);
  %absTol= quantile(uSlice(:),0.1); if absTol==0, absTol= defaultAbsTol; end
  absTol= params.viz.integrateAbstol;
  total_u(i)= fewcell.util.integrate(ufun,params.g.domainLim,absTol);
  i= i+1;
end
t= result.SolutionTimes(t);   % time idx -> real time
toc;
