function [u,x,y,t,total_u]= interpolateIntegrateAHL(model,result, ...
            params,startTimeIdx,stopTimeIdx)
          
% Precompute interpolation for all times and store it
if nargin<4, startTimeIdx= 1; stopTimeIdx= length(result.SolutionTimes); end
if nargin==4, error('[interpolateAHL]: specify end time'); end
resolution= params.viz.interpResolution;

% X,Y grid
if length(resolution)==1, resolution(2)= resolution(1); end
x= linspace(params.viz.domLim(1,1),params.viz.domLim(1,2), resolution(1));
y= linspace(-params.viz.domLim(2,1),-params.viz.domLim(2,2),resolution(2));
[X,Y]= meshgrid(x,y);

% Slice TimeDependentResults and interpolate + integrate each slice
tic;
timePoints= min(params.viz.timePoints,length(result.SolutionTimes));
t= floor(linspace(startTimeIdx,stopTimeIdx, timePoints));
u= zeros([size(X),length(t)]);
total_u= zeros([1,length(t)]);
if ~params.viz.figID(1), t= result.SolutionTimes(t); return; end
resultSlice= cell(timePoints,1);
for i= 1:length(t)
  resultSlice{i}= util.sliceResult(model,result,t(i));
end
intAbsTol= params.viz.integrateAbstol;
intDomainLim= params.g.domainLim;
parfor i= 1:length(t)
  ufun= @(x,y) util.nan0(reshape(interpolateSolution(resultSlice{i}, x,y), size(x)));
  u(:,:,i)= ufun(X,Y);
  % rectangular domain (r,z) & CYLINDRICAL coords
  total_u(i)= fewcell.util.integrate(ufun,intDomainLim,intAbsTol);
end
t= result.SolutionTimes(t);   % time idx -> real time
toc;
end
