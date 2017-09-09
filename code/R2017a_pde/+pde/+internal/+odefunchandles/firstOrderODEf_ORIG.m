function f=firstOrderODEf(t,u)
%firstOrderODEf - Residue function for first-order ODE system
%Computes residue of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2016 The MathWorks, Inc.

global femodel


if ~(femodel.vq || femodel.vg || femodel.vh || femodel.vr || ...
        femodel.vc || femodel.va || femodel.vf)
    f=-femodel.K*u+femodel.F;
    return
end

if(femodel.numConstrainedEqns == femodel.totalNumEqns)
    uFull = u;
else
  %!
    uFull = femodel.B*u + femodel.ud;
end

if femodel.vq || femodel.vg || femodel.vh || femodel.vr
  %!
    [femodel.Q,femodel.G,femodel.H,femodel.R]=femodel.thePde.assembleBoundary(uFull,t);
end

if femodel.vc || femodel.va || femodel.vf
    if(femodel.nrp==2)
        [femodel.K, femodel.F] = formGlobalKF2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t);
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
  %!
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
      %!
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

f=-K*u+F;

end
