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

% Singlecell Constants
SINormsM=1/(60e-9); SINorms=1/60; 
k1=0.01*SINormsM;
dAHL=0.01*SINorms;

global bactNodes;
N= length(bactNodes);
ahl= mean(u(bactNodes));
modelJacobian= singlecell.modelJacobian(1,[u(end-7:end-3);ahl;u(end-2:end)]);
% N: pde mesh nodes, y: singlecell variables
% 1: dn/dn
dndn= df;
dndn(bactNodes,:)= dndn(bactNodes,:) - repmat(mean(dndn(bactNodes,:),1),N,1);
%dndn(bactNodes,bactNodes)= dndn(bactNodes,bactNodes) - (k1*u(end-3)+dAHL)/N;
x= dndn(bactNodes,bactNodes);
A= eye(N)-1/N; A= A(1:end-1,:);
%Ainv= ones(N-1)+eye(N-1);
x0= x;
for i= 1:N
  x(:,i)= A\x0(1:end-1,i);
  %x(1:end-1,i)= Ainv*x0(1:end-1,i);
end
%x(end,:)= 0;%x0(end,:);
sum(sum(abs(x-x0-repmat(mean(x),N,1)))) ./N^2
svds(x,1,'smallest')
dndn(bactNodes,bactNodes)= x;

% 2: dn/dy
dndy= zeros(size(df,1),8);
dndy(bactNodes,:)= repmat(modelJacobian(6,[1:5,7:9]), N,1);
% 3: dy/dy
dydy= modelJacobian([1:5,7:9],[1:5,7:9]);
% 4: dy/dn
dydn= zeros(8, size(df,2));
dydn(:,bactNodes)= repmat(modelJacobian([1:5,7:9],6) ./ N, 1,N);
% Concatenate the pieces
df= [dndn, dndy; dydn, dydy];

%min(svd(dydy))
%det(dydy)
end
