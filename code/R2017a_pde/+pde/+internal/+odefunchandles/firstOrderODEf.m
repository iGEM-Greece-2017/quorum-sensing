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
global bactNodeIdx;
global bactProdMultiplier;  % multiplies each bacterium's AHL change, to simulate a denser bacterial population
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
    femodel.A = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,u(nodeIdx,i),t,'a');
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

[modelP,modelGrowth]= singlecell.modelCoeffs_weber(ones(11,1),1,false,false);
yIdx1= max(nodeIdx)+1;  % Index of 1st y of 0th bacterium
if bactNodesEqulength
  bactNodeN= sum(bactNodes(:,1));
  uCoeffNorm= F(bactNodeIdx);
  uCoeffNorm= uCoeffNorm ./ repmat( sum(uCoeffNorm) ,bactNodeN,1);
end

f= u;   % allocate only
for i= 1:size(u,2)
  yBact= zeros(11,nBact); yBact(11,:)= 1;
  yBact([1:5,7:9],:)= reshape(u(yIdx1: yIdx1+8*(nBact-1)+7, i), 8,nBact);
  % Calculate average AHL level for each bacterium
  if bactNodesEqulength
    uBact= reshape(u(bactNodeIdx,i), [],nBact);
    yBact(6,:)= sum(uBact.*uCoeffNorm);
  else
    % Stripped for efficiency in parfor
    %
    for b= 1:nBact
      uCoeff= F(bactNodes(:,b));
      yBact(6,b)= sum(u(bactNodes(:,b),i).*uCoeff./sum(uCoeff));    % y variables for this bacterium
    end
    %}
  end
  dy= singlecell.model_weber(t,yBact,modelP,modelGrowth);
  ahlProd= F;
  if bactNodesEqulength
    ahlProd(bactNodeIdx)= repmat(bactProdMultiplier.*dy(6,:), bactNodeN,1) .* F(bactNodeIdx);
  else
    %
    for b= 1:nBact
      ahlProd(bactNodes(:,b))= bactProdMultiplier(b)*dy(6,b)*F(bactNodes(:,b));
    end
    %}
  end
  f(:,i)= [ahlProd; reshape(dy([1:5,7:9],:),[],1)];
end
f= f + [-K*u(nodeIdx,:); zeros(nBact*8,size(u,2))];
end
