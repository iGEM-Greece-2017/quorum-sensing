function f=firstOrderODEf(t,u)
%firstOrderODEf - Residue function for first-order ODE system
%Computes residue of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2016 The MathWorks, Inc.

global femodel
%{ @@@ Custom @@@ %%
global enableSinglecellEq;
global bactNodes;
global bactWallDilution;
nBact= length(bactNodes);
%} @@@ Custom @@@ %%


if ~(femodel.vq || femodel.vg || femodel.vh || femodel.vr || ...
        femodel.vc || femodel.va || femodel.vf)
  if enableSinglecellEq, u= u(1:end-nBact*8,:); end
  f=-femodel.K*u+femodel.F*ones(1,size(u,2));
  return;
end

% make uFull
if(femodel.numConstrainedEqns == femodel.totalNumEqns)
  if enableSinglecellEq, uFull= u(1:end-nBact*8,:);
  else, uFull= u; end
else
  error('[firstOrderODEf]: code has been stripped!');
end

if femodel.vq || femodel.vg || femodel.vh || femodel.vr
  error('[firstOrderODEf]: code has been stripped!');
end

% [K,F,A]= formGlobalKF2D(uFull)
if femodel.vc || femodel.va || femodel.vf
  if(femodel.nrp==2)
          %{ ### @@@    DANGER, i==1 drops most of the computations !!!    @@@ ###
    i= 1;
    [femodel.K, femodel.F]= ...
      formGlobalKF2D(femodel.emptyPDEModel, femodel.p, ...
                     femodel.t, femodel.coefstruct,uFull(:,i),t);
    A = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull(:,i),t,'a');
    femodel.A = A; % TEMP - Change the femodel to have a A matrix later.
          %} ### @@@    DANGER, i==1 drops most of the computations !!!    @@@ ###
  elseif(femodel.nrp==3)
    error('[firstOrderODEf]: code has been stripped!');
  end
end

if ~(femodel.vh || femodel.vr)
  % neither H nor R are functions of t or u
  K = femodel.K + femodel.A + femodel.Q;
  F = femodel.B'*(femodel.F + femodel.G - K*femodel.ud);
  K = femodel.B'*K*femodel.B;
else
  error('[firstOrderODEf]: code has been stripped!');
end

f= zeros(size(u));
%% @@@ Custom @@@ %%
if ~enableSinglecellEq    
  error('[firstOrderODEf]: code has been stripped!');
end
% Initial calculation for d(n_i)/dt
f= [-K*uFull; zeros(nBact*8,size(u,2))];

  %global dyLog;

for i= 1:size(u,2)
  for b= 1:nBact
    physToComp= F(bactNodes{b});
    ahl= sum(u(bactNodes{b},i).*physToComp./sum(physToComp));
    
    yIdx= (nBact-b+1)*8-1;    % where the y variables for this bacterium start (relative to end)
    yBact= [u(end-yIdx: end-yIdx+4,i); ahl; u(end-yIdx+5: end-yIdx+7,i); 0;1];   % y variables for this bacterium

    dy= singlecell.model_weber(t,yBact,false,false);
    ahlProd= dy(6)*physToComp;

    % dy(6) contains the total derivative of all the contributing terms
    % dy(6):= d(ahl)/dt
    % From the definition of <ahl>, an equation that must hold can be derived with
    % the chain rule: mean(n_i')==AHL'
    % A correction factor will be added to the n_i' to accomodate it:

    % Note:
    %   The cylindrical symmetry allows for many bacteria that sit on a cylindrical ring
    % of set radius to be simulated by a single set of singlecell equations. Their total
    % effect is incorporated by multiplying their number with their output (dy_6)
    f(bactNodes{b},i)= f(bactNodes{b},i) + ahlProd;

    f(end-yIdx: end-yIdx+7,i)= dy([1:5,7:9]);
  end

  %log10([max(u(1:end-16)) dy(6)])

  % Don't allow nodes to stay negative
  %f(u<=0 & f<=0)= -f(u<=0 & f<=0);%-u(u<=0 & f<=0);
end
end
