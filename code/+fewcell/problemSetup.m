function [p,model,tlist,domainVolume,bactRingDensity]= problemSetup(p,plotMesh)   %param
% Setup all elements of the pde problem

%% Time
%ntime= (p.t.tstop - p.t.tstart)+1;
tlist= linspace(p.t.tstart, p.t.tstop, p.t.timePoints);
%% Geometry
domainVolume= pi*p.g.domainLim(1).^2 * p.g.domainLim(2);  % domain is actually a cylinder
bactVolume= pi/4*p.g.bactSize(1).^2 * p.g.bactSize(2);
% Develop bactCenters from metaparameters
p.g.bactCenters= zeros(p.g.nRings*p.g.nLayers,2);
p.g.bactCenters(1,:)= p.g.bactCenter0;
prevLayerStart= 1;
for layer= 1:p.g.nLayers
  if layer>1
    p.g.bactCenters((layer-1)*p.g.nRings+1,:)= ...
      p.g.bactCenters(prevLayerStart,:)-[0,p.g.layerSeparation+p.g.bactSize(2)];
  end
  for ring= (layer-1)*p.g.nRings+2:layer*p.g.nRings
    p.g.bactCenters(ring,:)= ...
      p.g.bactCenters(ring-1,:)+[p.g.ringDist+p.g.bactSize(1),0];
  end
  prevLayerStart= (layer-1)*p.g.nRings+1;
end
nBact= size(p.g.bactCenters,1);
global bactProdMultiplier; bactProdMultiplier= p.g.bactProdMultiplier;

bactRingDensity= fewcell.util.bactRingDensity(p.g.bactCenters(1:p.g.nRings,1),p.g.bactSize, 0);
bactRingDensity= bactRingDensity*p.g.nLayers*bactProdMultiplier;
totalBacteria= sum(bactRingDensity);
if p.growth.on && p.growth.maxRings < p.g.nRings
  totalBacteria= sum(bactRingDensity(end-p.growth.maxRings:end));
end
fprintf('Max bacteria: %d', round(totalBacteria));

% Create geometry description
x= [1,0;0,1;0,1;1,0]; y= [1,0;1,0;0,1;0,1];
names= cell(2,nBact+1);
names(:,1)= {'domain', '+'};
shapes= zeros(10,nBact+1);
shapes(:,1)= [3;4; x*[0;p.g.domainLim(1)]; y*[0;-p.g.domainLim(2)]];  %domain
% Bacterial geometry
for b= 1:nBact
  bsize= p.g.bactSize;
  bULcorner= p.g.bactCenters(b,:) - [bsize(1)/2,-bsize(2)/2];   % upper left corner
  shapes(:,b+1)= [3;4; x*[bULcorner(1); bULcorner(1)+bsize(1)]; y*[bULcorner(2); bULcorner(2)-bsize(2)]];
  names(:,b+1)= {['bact',num2str(b)]; '+'};
end
% Create the description
names{2,end}= '';
setf= [names{:}];
names= char(names{1,:}); names= names';
[geometryDescription,~]= decsg(shapes,setf,names);

%% Growth
if p.growth.on, fprintf('\tGrowth ON\n'); else, fprintf('\tGrowth OFF\n'); end
% Calculate growth step sizes
global growth;
global enableSinglecellEq;
global enableGraphics
if ~enableSinglecellEq, p.growth.on= false; end
growth.on= p.growth.on;
if growth.on
  % Growth curve params
  p.growth.params.Nmax= totalBacteria;
  p.growth.params.Nmin= 0;   % disregard lag phase
  % Copy parameters to global
  growth.r0= p.growth.r0; growth.dr= p.growth.dr;
  growth.maxRings= p.growth.maxRings; growth.nLayers= p.g.nLayers;
  growth.initDNA= p.solve.y0(1); growth.bactRingDensity= bactRingDensity;
  maxStep= floor(p.growth.maxRings/p.growth.dr);
  
  % Sum the number of bacteria every <dr> rings (for all layers) to calculate the growth step size
  assert(growth.dr<=growth.r0, 'Number of rings to grow at each growth step must be <= that the initial rings');
  assert(~mod(p.g.nRings-growth.r0, growth.dr), '<nRings>-<r0> is not divisible by <dr>; increase <nRings>');
  growth.stepSize= sum(reshape(bactRingDensity(growth.r0+1:end),p.growth.dr,[]),1);
  if p.growth.maxRings < p.g.nRings
    growth.stepSize= growth.stepSize - [zeros(1,maxStep),growth.stepSize(1:end-maxStep)];
  end
  % Calculate growth curve
  growth.bactN= singlecell.growthCurve(sum(bactRingDensity(1:p.growth.r0)),tlist,p.growth.params);
  % Calculate the time at which each growth step should occur
  [growth.tstep, tLast]= fewcell.util.makeGrowthTimesteps(growth);

  % Plot growth curve vs quantized growth curve
  qgc= zeros(size(growth.bactN))+sum(reshape(bactRingDensity(1:growth.r0),p.growth.dr,[]),1);
  for t= 1:length(growth.tstep)-1
    qgc(growth.tstep(t):growth.tstep(t+1))= qgc(growth.tstep(t))+growth.stepSize(t);
  end
  if plotMesh && enableGraphics
    figure(1); plot(tlist/60,[growth.bactN,qgc]);
  end
else
  growth.tstep= length(tlist);
end
growth.minTstep= min(tlist(growth.tstep) - tlist([1 growth.tstep(1:end-1)]));
if growth.on, fprintf('Num of growth timesteps: %d\tMin timestep: %f\n', length(growth.tstep),growth.minTstep); end

%% PDE
model= createpde(1);
geometryFromEdges(model,geometryDescription);
%% Boundaries
numEdges= model.Geometry.NumEdges;
applyBoundaryCondition(model, 'neumann', 'Edge',1:numEdges,'g',0,'q',0);
%% Coefficients
% The del operator is expressed in cylindrical coordinates. du/dθ=0 due to symmetry, so:
%    u' - del*( c del(u)) +  au =  f
% becomes, after substituting r<-x,
%   xu' - del*(xc del(u)) + xau = xf
% c: diffusion, a: destruction rate
dCoeff= @(r,s)r.x;
cCoeff= @(r,s)p.c.c_agar*r.x;
aCoeff= @(r,s)p.c.d_AHL*r.x;
specifyCoefficients(model,'Face',1,'m',0,'d',dCoeff,'c',cCoeff,'a',aCoeff,'f',0);
for b=1:nBact
  %bactProd= @(r,s)(10*s.u+(1e1*(s.time<2)));
  %bactProd= @(r,s)(1*(s.time==tstart));
  %bactProd= @(r,s)1e4*(r.x+s.u);
  % Spatial coefficient for production through singlecell eqs
  %bactProd= @(r,s)2*pi*r.x*p.g.bactSize(1)*p.g.bactSize(2)./bactVolume;
  bactProd= @(r,s)r.x;
  cCoeff= @(r,s)p.c.c_cytoplasm*r.x;
  specifyCoefficients(model,'Face',b+1,...
    'm',0,'d',dCoeff,'c',cCoeff,'a',0,'f',bactProd);
end
%% Initial conditions
setInitialConditions(model,0);
%% Mesh
generateMesh(model,'MesherVersion','R2013a','GeometricOrder','linear', 'Jiggle','minimum','JiggleIter',50,...
  'Hgrad',p.m.Hgrad, 'Hmax',p.g.domainLim(1) * p.m.HmaxCoeff);
totalMeshNodes= size(model.Mesh.Nodes,2);
fprintf('Total mesh nodes: %d\n', totalMeshNodes);
fprintf('Total equations: %d\n', totalMeshNodes + nBact*8);

% Plot mesh
if plotMesh && enableGraphics
  figure(2); clf;
  pdeplot(model,'NodeLabels','on');
  %pdegplot(model, 'EdgeLabels','on','FaceLabels','on');
  axis tight;
  meshPlot= gca;
  meshPlot.XLim= p.viz.domLim(1,:);
  meshPlot.YLim= -flip(p.viz.domLim(2,:));
  drawnow;
end
