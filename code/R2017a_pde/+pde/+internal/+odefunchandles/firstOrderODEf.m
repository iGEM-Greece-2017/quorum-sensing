function f=firstOrderODEf(t,u)
%firstOrderODEf - Residue function for first-order ODE system
%Computes residue of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2016 The MathWorks, Inc.

global femodel
%{ @@@ Custom @@@ %%
global enableSinglecellEq;
global bactNodesEqulength;
global bactNodes;
%nBact= length(bactNodes);
nBact= size(bactNodes,2);
%} @@@ Custom @@@ %%


if ~(femodel.vq || femodel.vg || femodel.vh || femodel.vr || ...
        femodel.vc || femodel.va || femodel.vf)
  if enableSinglecellEq, u= u(1:end-nBact*8,:); end
  f=-femodel.K*u+femodel.F*ones(1,size(u,2));
  return;
end

% make uFull
nodeIdx= 1:size(u,1)-enableSinglecellEq*nBact*8;
if(femodel.numConstrainedEqns == femodel.totalNumEqns)
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
                     femodel.t, femodel.coefstruct,u(nodeIdx,i),t);
    A = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,u(nodeIdx,i),t,'a');
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

%% @@@ Custom @@@ %%
if ~enableSinglecellEq    
  error('[firstOrderODEf]: code has been stripped!');
end
% Initial calculation for d(n_i)/dt
f= [-K*u(nodeIdx,:); zeros(nBact*8,size(u,2))];

[modelP,modelGrowth]= singlecell.modelCoeffs_weber(ones(11,1),false,false);
yIdx1= max(nodeIdx)+1;  % Index of 1st y of 0th bacterium


for i= 1:size(u,2)
  yBact= zeros(11,nBact); yBact(11,:)= 1;
  yBact([1:5,7:9],:)= reshape(u(yIdx1: yIdx1+8*(nBact-1)+7, i), 8,nBact);
  for b= 1:nBact
    physToComp= F(bactNodes(:,b));
    yBact(6,b)= sum(u(bactNodes(:,b),i).*physToComp./sum(physToComp));    % y variables for this bacterium
  end
  dy= singlecell.model_weber(t,yBact,modelP,modelGrowth);
  ahlProd= F;
  for b= 1:nBact
    ahlProd(bactNodes(:,b))= dy(6,b)*F(bactNodes(:,b));
    %physToComp= F(bactNodes(:,b));
    %ahlProd= dy(6,b)*physToComp;

    % dy(6) contains the total derivative of all the contributing terms
    % dy(6):= d(ahl)/dt
    % From the definition of <ahl>, an equation that must hold can be derived with
    % the chain rule: mean(n_i')==AHL'
    % A correction factor will be added to the n_i' to accomodate it:

    % Note:
    %   The cylindrical symmetry allows for many bacteria that sit on a cylindrical ring
    % of set radius to be simulated by a single set of singlecell equations. Their total
    % effect is incorporated by multiplying their number with their output (dy_6)
    %f(bactNodes(:,b),i)= f(bactNodes(:,b),i) + ahlProd;
  end
  f(1:yIdx1-1,i)= f(1:yIdx1-1,i) + ahlProd;
  f(yIdx1: yIdx1+8*(nBact-1)+7, i)= reshape(dy([1:5,7:9],:),[],1);
end
end
