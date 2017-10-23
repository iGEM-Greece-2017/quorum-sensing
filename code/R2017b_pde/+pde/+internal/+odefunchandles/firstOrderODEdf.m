function df=firstOrderODEdf(t,u)
%firstOrderODEdf - Jacobian function for first-order ODE system
%Computes jacobian of discretized PDE with first-order time derivative.
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

if ~(femodel.vq || femodel.vg || femodel.vh || femodel.vr || femodel.vc || femodel.va || femodel.vf)
  df=-femodel.K;
  return
end

% Make uFull
nodeIdx= 1:size(u,1)-enableSinglecellEq*nBact*8;
if(femodel.numConstrainedEqns == femodel.totalNumEqns)
else
  error('[firstOrderODEdf]: code has been stripped!');
end

if(femodel.vq || femodel.vh)
  error('[firstOrderODEdf]: code has been stripped!');
end

if(femodel.vc || femodel.va)
  if(femodel.nrp==2)
    i= 1;
    [K, ~]= formGlobalKF2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,u(nodeIdx,i),t);
    A= formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,u(nodeIdx,i),t,'a');
  elseif(femodel.nrp==3)
    error('[firstOrderODEdf]: code has been stripped!');
  end
end
if femodel.vq, femodel.Q=Q; end
if femodel.vh
  error('[firstOrderODEdf]: code has been stripped!');
end
if femodel.vc, femodel.K=K; end
if femodel.va, femodel.A=A; end
if ~femodel.vh
  df=-femodel.B'*(femodel.K+femodel.A+femodel.Q)*femodel.B;
else
  error('[firstOrderODEdf]: code has been stripped!');
end

% N: pde mesh nodes, y: singlecell variables
% 1: dn'/dn
% -> Remains unchanged (df)
% 2: dn'/dy
dndy= zeros(length(nodeIdx),8*nBact);
% 3: dy'/dn
dydn= zeros(8*nBact,length(nodeIdx));
% 4: dy'/dy
dydy= zeros(8*nBact);

[modelP,modelGrowth]= singlecell.modelCoeffs_weber(ones(11,1),false,false);
yIdx1= max(nodeIdx)+1;  % Index of 1st y of 0th bacterium
if bactNodesEqulength
  bactNodeN= sum(bactNodes(:,1));
  uCoeffNorm= femodel.F(bactNodeIdx);
  uCoeffNorm= uCoeffNorm ./ repmat( sum(uCoeffNorm) ,bactNodeN,1);
end
for i= 1:1%size(u,2)
  yBact= zeros(11,nBact); yBact(11,:)= 1;
  yBact([1:5,7:9],:)= reshape(u(yIdx1: yIdx1+8*(nBact-1)+7, i), 8,nBact);
  % Calculate average AHL level for each bacterium
  if bactNodesEqulength
    uBact= reshape(u(bactNodeIdx,i), [],nBact);
    yBact(6,:)= sum(uBact.*uCoeffNorm);
  else
    for b= 1:nBact
      uCoeff= femodel.F(bactNodes(:,b));
      yBact(6,b)= sum(u(bactNodes(:,b),i).*uCoeff./sum(uCoeff));    % y variables for this bacterium
    end
  end
  j= singlecell.modelJacobian_weber(t,yBact,modelP,modelGrowth);
  for b= 1:nBact
    bIdx= (b-1)*8+1:b*8;
    dndy(bactNodes(:,b),bIdx)= femodel.F(bactNodes(:,b)).*bactProdMultiplier(b).*j(6,[1:5,7:9],b);
    dydn(bIdx,bactNodes(:,b))= j([1:5,7:9],6,b)*femodel.F(bactNodes(:,b))';
    dydy(bIdx, bIdx)= j([1:5,7:9],[1:5,7:9],b);
  end
end

df= [df,dndy; dydn,dydy];
%min(svd(dydy))
%det(dydy)
%spy(df);
end
