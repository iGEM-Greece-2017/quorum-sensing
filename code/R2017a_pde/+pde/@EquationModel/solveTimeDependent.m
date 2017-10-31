function [u,dudt] = solveTimeDependent(self, coefstruct, u0, ut0, tlist, tsecondOrder)
% solveTimeDependent - Internal function for solving a time-dependent PDE
%

% Copyright 2015-2016 The MathWorks, Inc.

global femodel
clearN = onCleanup(@() clear('global', 'femodel'));
import pde.internal.odefunchandles.*

rtol =  self.SolverOptions.RelativeTolerance;
atol =  self.SolverOptions.AbsoluteTolerance;
stats = self.SolverOptions.ReportStatistics;
if ~tsecondOrder
    dudt = [];
end

[p,e,t] = self.Mesh.meshToPet();
femodel=pde.DynamicDiscretizedPDEModel(self,p,e,t,coefstruct,u0,tlist,tsecondOrder);
%{ @@@ Custom @@@ %%
  global enableSinglecellEq;
  global solveInternalParams;
  global bactNodes;
  global yyResults;
  
  if enableSinglecellEq, femodel.vf= 1; end
  %nBact= length(bactNodes);
  nBact= size(bactNodes,2);
%} @@@ Custom @@@ %%

if(~femodel.IsSpatialCoefficientsNonlinear)
    % Spatial coefficients are independent of solution and time.
    % Impose BC on the system matrices only once in the beginning.
    nun=femodel.numConstrainedEqns;
    [null,orth]=pdenullorth(femodel.H);
    if size(orth,2)==0
        ud=zeros(size(femodel.K,2),1);
    else
        ud=full(orth*((femodel.H*orth)\femodel.R));
    end
    KK=femodel.K+femodel.A+femodel.Q;
    FF=null'*((femodel.F+femodel.G)-femodel.K*ud);
    KK=null'*KK*null;
    KK = self.checkForSymmetry(KK);
    if femodel.IsSecondOrderODE
        femodel.K=[sparse(nun,nun) -speye(nun,nun); KK sparse(nun,nun)];
        femodel.F=[zeros(nun,1);FF];
    else
        femodel.K=KK;
        femodel.F=FF;
    end
end

nu=size(p,2)*self.PDESystemSize;
nuu=nu;

if tsecondOrder
    [u0,ut0] = imposeBConICsecondOrderODE(self,u0,ut0,tlist);
    nu=size(u0,1);
    uu0=[u0;ut0];
    if (self.EquationCoefficients.dDefined())
        femodel.C = coefstruct.Coefficients.d{1};
    end
    %Following pragma, %#function, are needed for deployment.
    %#function pde.internal.odefunchandles.secondOrderODEf
    fcn  = str2func('secondOrderODEf');
    %#function pde.internal.odefunchandles.secondOrderODEdf
    dfcn = str2func('secondOrderODEdf');
    %#function pde.internal.odefunchandles.secondOrderODEm
    mfcn = str2func('secondOrderODEm');
    odeoptions = constructODEoptions(uu0,rtol,atol,stats,nu,dfcn,mfcn,tsecondOrder);
else
    u0 = imposeBConICfirstOrderODE(self,u0,tlist);
    nu=size(u0,1);
    uu0 = u0;
    %Following pragma, %#function, are needed for deployment.
    %#function pde.internal.odefunchandles.firstOrderODEf
    fcn  = str2func('firstOrderODEf');
    %#function pde.internal.odefunchandles.firstOrderODEdf
    dfcn = str2func('firstOrderODEdf');
    %#function pde.internal.odefunchandles.firstOrderODEm
    mfcn = str2func('firstOrderODEm');
    odeoptions = constructODEoptions(u0,rtol,atol,stats,nu,dfcn,mfcn,tsecondOrder);
end

%{ @@@ Custom @@@ %%
  
  nNode= length(u0);
  % ODE solver options
  odeoptions.Vectorized= 'on';
  odeoptions.Jacobian= [];
  odeoptions.JConstant= 'off';
  odeoptions.InitialStep= 2;
  if enableSinglecellEq
    nY= nBact*8;  % number of extra variables
    odeoptions.AbsTol= [odeoptions.AbsTol*ones(1,length(uu0)), zeros(1,nY)];
    odeoptions.AbsTol(end-(nY-1):end)= repmat(solveInternalParams.AbsTol_y, 1,nBact);
    uu0= [uu0; solveInternalParams.y0];
  end
  
  %% Sparse Jacobian
  %testInput= 2*rand(size(uu0))+1;
  testInput= ones(size(uu0));
  JPattern= (speye(length(uu0))==1) | pde.internal.odefunchandles.firstOrderODEdf(1,testInput);   % init
  %
  global bactNodesEqulength; global bactNodeIdx; global bactProdMultiplier;
  fcnJ= @(t,u)pde.firstOrderODEf_noGlobal(t,u,femodel,true,bactNodesEqulength,bactNodes,bactNodeIdx,bactProdMultiplier);
  %dfJ= @(t,u)pde.firstOrderODEdf_noGlobal(t,u,femodel,enableSinglecellEq,bactNodesEqulength,bactNodes,bactNodeIdx,bactProdMultiplier);
  jOpt= struct('diffvar',2,'vectvars',2,'thresh',odeoptions.AbsTol','fac',[]);
  %nElt= [0,10*16+1];
  %while nElt(end) - nElt(1) > 10*16
    parfor i= 1:2
      if i==1         % QS transition start
        test_y0= [1.365;4.186;1.137;196.5;621.7; 23.18;25.00;0.1697];
        test_AHL= 3.801;
      else            % QS transition end
        test_y0= [0.0281;43.40;10.85;5897;18800; 466.0;10740;1.5066];
        test_AHL= 2.603;
      end
      testInput= [test_AHL*ones(size(u0))+rand(size(u0)); repmat(test_y0,nBact,1)]; 
      JPattern= JPattern | (pde.odenumjac(fcnJ,{1,testInput},fcnJ(1,testInput),jOpt) ~= 0);
      %JPattern= JPattern | dfJ(1,testInput);
    end
    %JPattern(1:nNode,1:nNode)= JPatternNode;
    %nElt= [nElt(2:end),nnz(JPattern)]
  %end
  %}
  odeoptions.JPattern= JPattern;
  %% Solve
  solution= ode15s(fcn,tlist,uu0,odeoptions);
  assert(tlist(end) <= solution.x(end));
  uu= deval(solution,tlist)';
  if enableSinglecellEq
    yyResults= cell(nBact,1);
    for b= 1:nBact
      % where the y variables for this bacterium start (relative to end)
      yIdx= (nBact-b+1)*8-1;
      ahl= mean(uu(:,bactNodes(:,b)),2);
      % y variables for this bacterium
      yyResults{b}= [uu(:,end-yIdx: end-yIdx+4), ahl, uu(:,end-yIdx+5: end-yIdx+7), ahl];
    end
    uu= uu(:,1:end-nY);
  end
  
%} @@@ Custom @@@ %%

numCompTimes = size(uu, 1);
if(numCompTimes < size(tlist, 2))
    warning(message('pde:pdetool:incompSoln'))
end

DoFIndex=1:nu; % vector of 1 to number of degree of freedom
if(length(tlist)==2)
    % ODE15S behaves differently if length(tlist)==2
    numCompTimes = 2;
    u1(DoFIndex,:)    = transpose(uu([1 size(uu,1)],DoFIndex));
else
    u1(DoFIndex,:)    = transpose(uu(:,DoFIndex));
end

du1dt = zeros(size(u1));

if tsecondOrder
    velIndex=nu+1:size(uu,2);
    if(length(tlist)==2)
        % ODE15S behaves differently if length(tlist)==2
        numCompTimes = 2;
        du1dt(DoFIndex,:) = transpose(uu([1 size(uu,1)],velIndex));
    else
        du1dt(DoFIndex,:) = transpose(uu(:,velIndex));
    end
end
%
% If Dirichlet constraints are functions of time or u, recovery of the
% full solution requires reevaluating H and R matrices
if (femodel.vh || femodel.vr)
    % new constrained values calc
    uu1=zeros(nuu,numCompTimes);
    duu1=zeros(nuu,numCompTimes);
    for k=1:numCompTimes
        [uu1(:,k), duu1(:,k)] = recoverFullSolution(self,u1(:,k),du1dt(:,k),tlist(k),femodel.tspan,femodel.H,femodel.R, ...
            rtol, atol(1));
    end
    u1=uu1;
    if tsecondOrder
        du1dt = duu1;
    end
else
    u1=femodel.B*u1+femodel.ud*ones(1,numCompTimes);
    if tsecondOrder
        du1dt=femodel.B*du1dt+femodel.dudt*ones(1,numCompTimes);
    end
end
u = u1;
if tsecondOrder
    dudt = du1dt;
end
end


function u0 = imposeBConICfirstOrderODE(thePde, u0,tlist)
global femodel
if ~(femodel.vh || femodel.vr)
    u0=femodel.B'*u0;
else
    [~,~,H,~]=thePde.assembleBoundary(u0,tlist(1));
    B=pdenullorth(H);
    u0=B'*u0;
end

end

function [u0, ut0] = imposeBConICsecondOrderODE(thePde, u0,ut0,tlist)
global femodel
if ~(femodel.vh || femodel.vr)
    u0=femodel.B'*u0;
    ut0=femodel.B'*ut0;
else
    t0=tlist(1);
    t1=t0+femodel.tspan*sqrt(eps);
    dt=t1-t0;
    [~,~,H,R]=thePde.assembleBoundary(u0,t0);
    [N,O]=pdenullorth(H);
    ud=O*((H*O)\R);
    [~,~,H,R]=thePde.assembleBoundary(u0,t1);
    [N1,O]=pdenullorth(H);
    ud1=O*((H*O)\R);
    u0=N'*u0;
    ut0=N'*(ut0-((N1-N)*u0-(ud1-ud))/dt);
end
end


function odeoptions = constructODEoptions(uu0,rtol,atol,stats,nu,dfcn,mfcn,tsecondOrder)
global femodel

odeoptions = odeset;
% If doJacFiniteDiff is true, ode15s computes the Jacobian by finite
% difference. This is useful for testing purposes.
doJacFiniteDiff=false;
if(doJacFiniteDiff)
    jac0=fcn(0,uu0);
    [ir,ic]=find(jac0);
    jpat=sparse(ir,ic,ones(size(ir,1),1), size(jac0,1), size(jac0,2));
    odeoptions=odeset(odeoptions, 'JPattern', jpat);
else
    odeoptions=odeset(odeoptions, 'Jacobian',dfcn);
end
if(~(femodel.vc || femodel.va || femodel.vq || femodel.vh || femodel.vr))
    odeoptions=odeset(odeoptions,'JConstant','on');
end

if(~(femodel.vd || femodel.vh))
    odeoptions=odeset(odeoptions,'Mass', mfcn(0, uu0));
else
    odeoptions=odeset(odeoptions,'Mass',mfcn);
end
odeoptions=odeset(odeoptions,'RelTol',rtol);

if (tsecondOrder)
    % Mask out dudt
    atol=atol*ones(2*nu,1);
    atol(nu+1:nu+nu)=Inf*ones(nu,1);
    odeoptions=odeset(odeoptions,'MaxOrder',2);
end

odeoptions=odeset(odeoptions,'AbsTol',atol);
odeoptions=odeset(odeoptions,'Stats',stats);

end


function [uf, dudtf, isConverged] = recoverFullSolution(thePde,ur,durdt,t,tspan,H,R, rtol, atol)
%This function is used to transform the reduced solution with constrained DOFs
% removed to the full solution with all DOFs at all FE nodes.
%
% When the H or R matrices are functions of the dependent variables
% (i.e. a nonlinear system), multiple iterations may be required to
% obtain a converged solution.

maxIter = 50;
ufi = zeros(size(H,2),1);

dt=tspan*eps^(1/3);
t1 = t+dt;
tm1 = t-dt;

for i=1:maxIter
    [null,orth]=pdenullorth(H);
    ud=full(orth*((H*orth)\R));
    uf=null*ur + ud;
    % Compute time-derivative using a central difference approximation.
    [~,~,H,R]=thePde.assembleBoundary(uf,t1);
    ud1=full(orth*((H*orth)\R));
    [~,~,H,R]=thePde.assembleBoundary(uf,tm1);
    udm1=full(orth*((H*orth)\R));
    udDot = (ud1 - udm1)./(2*dt);
    dudtf = null*durdt + udDot;
    
    err=max(abs(ufi-uf)); % Check convergence only on displacements
    isConverged = all(abs(ufi-uf) < max(rtol*abs(uf), atol));
    if(isConverged)
        break;
    end
    ufi=uf;
    [~,~,H,R]=thePde.assembleBoundary(uf,t);
end
if(~isConverged)
    warning(message('pde:pdeCalcFullU:unconverged', ...
        sprintf('%12.3e', err)));
end
end


