function f=firstOrderODEf(t,u)
%firstOrderODEf - Residue function for first-order ODE system
%Computes residue of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2016 The MathWorks, Inc.

global femodel
%{ @@@ Custom @@@ %%
global bactNodes;
nBact= length(bactNodes);
%} @@@ Custom @@@ %%


if ~(femodel.vq || femodel.vg || femodel.vh || femodel.vr || ...
        femodel.vc || femodel.va || femodel.vf)
    f=-femodel.K*u+femodel.F;
    return
end

if(femodel.numConstrainedEqns == femodel.totalNumEqns)
    uFull= u(1:end-nBact*8);
else
    uFull = femodel.B*u + femodel.ud;
end

if femodel.vq || femodel.vg || femodel.vh || femodel.vr
    [femodel.Q,femodel.G,femodel.H,femodel.R]=femodel.thePde.assembleBoundary(uFull,t);
end

if femodel.vc || femodel.va || femodel.vf
    if(femodel.nrp==2)
        [femodel.K, femodel.F]= ...
          formGlobalKF2D(femodel.emptyPDEModel, femodel.p, ...
                         femodel.t, femodel.coefstruct,uFull,t);
        A = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'a');
        femodel.A = A; % TEMP - Change the femodel to have a A matrix later.
    elseif(femodel.nrp==3)
        [femodel.K, femodel.F] = formGlobalKF3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t);
        A = formGlobalM3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'a');
        femodel.A = A;
    end
end

if ~(femodel.vh || femodel.vr)
    % neither H nor R are functions of t or u
    K = femodel.K + femodel.A + femodel.Q;
    F = femodel.B'*(femodel.F + femodel.G - K*femodel.ud);
    K = femodel.B'*K*femodel.B;
else
    dt=femodel.tspan*sqrt(eps);
    t1=t+dt;
    dt=t1-t;
    if femodel.vd
        if(femodel.nrp==2)
            femodel.Mass = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'d');
        elseif(femodel.nrp==3)
            femodel.Mass = formGlobalM3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'d');
        end
    end
    if femodel.vh
        [N,O]=pdenullorth(femodel.H);
        ud=O*((femodel.H*O)\femodel.R);
        [~,~,H,R]=femodel.thePde.assembleBoundary(uFull,t1);
        [N1,O1]=pdenullorth(H);
        ud1=O1*((H*O1)\R);
        K=N'*((femodel.K+femodel.A+femodel.Q)*N + femodel.Mass*(N1-N)/dt);
        F=N'*((femodel.F+femodel.G)-(femodel.K+femodel.A+femodel.Q)*ud - ...
            femodel.Mass*(ud1-ud)/dt);
    else
        HH=femodel.H*femodel.Or;
        ud=femodel.Or*(HH\femodel.R);
        [~,~,~,R]=femodel.thePde.assembleBoundary(uFull,t1);
        ud1=femodel.Or*(HH\R);
        K=femodel.Nu'*(femodel.K+femodel.A+femodel.Q)*femodel.Nu;
        F=femodel.Nu'*((femodel.F+femodel.G)-(femodel.K+femodel.A+femodel.Q)*ud - ...
            femodel.Mass*(ud1-ud)/dt);
    end
    femodel.ud=ud;
end

%% @@@ Custom @@@ %%
% Initial calculation for d(n_i)/dt
f= [-K*u(1:end-nBact*8); zeros(nBact*8,1)];
for b= 1:nBact
  ahl= mean(u(bactNodes{b}));
  yIdx= (nBact-b+1)*8-1;    % where the y variables for this bacterium start (relative to end)
  yBact= [u(end-yIdx: end-yIdx+4); ahl; u(end-yIdx+5: end-yIdx+7)];   % y variables for this bacterium
  model= singlecell.model(1,yBact);
  % model(6) contains the total derivative of all the contributing terms
  % model(6):= d(ahl)/dt
  % From the definition of <ahl>, an equation that must hold can be derived with
  % the chain rule: mean(n_i')==AHL'
  % A correction factor will be added to the n_i' to accomodate it:
  % corr_ahl= AHL' - mean(n_i')
  %f(bactNodes{b}(1:2))= f(bactNodes{b}(1:2)) - mean(f(bactNodes{b}(1:2))) + model(6);
  f(bactNodes{b})= f(bactNodes{b}) - mean(f(bactNodes{b})) + model(6);

  f(end-yIdx: end-yIdx+7)= model([1:5,7:9]);
end

%{
% Notes for the Jacobian
  ni'= k-mean(k) + ahl'
% 1: dn'/dn
dni'/dnj= dk/dn + d(corr_ahl)/dn= dk/dn + d(ahl')/dn - 1/N*d(Σki)/dnj
= ki_nj + 1/N*(sum_i(dni'/dnj) - sum_i(ki_nj))

% 2: dy/dn
d(yi')/dnj= d(yi')/dAHL dAHL/dnj= d(yi')/dAHL * 1/N

% 3: dn/dy
dni'/dyj= dk/dyj - 1/N*sum_i(ki_yj) + d(ahl')/dyj - 
%}
end
