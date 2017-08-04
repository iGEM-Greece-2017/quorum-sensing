function df=firstOrderODEdf(t,u)
%firstOrderODEdf - Jacobian function for first-order ODE system
%Computes jacobian of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2016 The MathWorks, Inc.

global femodel

if ~(femodel.vq || femodel.vg || femodel.vh || femodel.vr || femodel.vc || femodel.va || femodel.vf)
  df=-femodel.K;
  return
end

if(femodel.numConstrainedEqns == femodel.totalNumEqns)
  uFull = u;
else
  uFull = femodel.B*u + femodel.ud;
end

if(femodel.vq || femodel.vh)
    [Q,~,H,~]=femodel.thePde.assembleBoundary(uFull,t);
end

if(femodel.vc || femodel.va)
    if(femodel.nrp==2)
        [K, ~] = formGlobalKF2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t);
        A = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'a');
    elseif(femodel.nrp==3)
        [K, ~] = formGlobalKF3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t);
        A = formGlobalM3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'a');
    end
end

if femodel.vq
  femodel.Q=Q;
end

if femodel.vh
  femodel.B=pdenullorth(H);
end

if femodel.vc
  femodel.K=K;
end

if femodel.va
  femodel.A=A;
end

if ~femodel.vh
  df=-femodel.B'*(femodel.K+femodel.A+femodel.Q)*femodel.B;
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
  [~,~,H,~]=femodel.thePde.assembleBoundary(uFull,t1);
  N1=pdenullorth(H);
  df=-femodel.B'*((femodel.K+femodel.A+femodel.Q)*femodel.B + ...
    femodel.Mass*(N1-femodel.B)/dt);
end


global bactNodes;
N= length(bactNodes);
ahl= mean(u(bactNodes));
modelJacobian= singlecell.modelJacobian(1,[u(end-7:end-3);ahl;u(end-2:end)]);
% N: pde mesh nodes, y: singlecell variables
% 1: dn/dn
dndn= df;
dndn(bactNodes(1),:)= 0;
% 2: dn/dy
dndy= zeros(size(df,1),8);
dndy(bactNodes(1),:)= modelJacobian(6,[1:5,7:9]);
% 3: dy/dy
dydy= modelJacobian([1:5,7:9],[1:5,7:9]);
% 4: dy/dn
dydn= zeros(8, size(df,2));
dydn(:,bactNodes(1))= modelJacobian([1:5,7:9],6);
% Concatenate the pieces
df= [dndn, dndy; dydn, dydy];

%min(svd(dydy))
%det(dydy)
spy(df);
end
