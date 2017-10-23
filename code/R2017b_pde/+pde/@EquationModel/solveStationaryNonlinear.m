function u = solveStationaryNonlinear(self, coefstruct,u0)
% solveStationaryNonlinear - Internal function for solving nonlinear stationary PDE
%

% Copyright 2015-2017 The MathWorks, Inc.

[p,~,t] = self.Mesh.meshToPet();
N = self.PDESystemSize;
[ndim, np] = size(p);

[K, A, F, Q, G, H, R] = assemblePdeMatrices(self, ndim, p, t, coefstruct, u0);

if ~any(any(H))
    H=sparse(0,size(H,2));
    R=sparse(0,1);
end

if(all(u0==0))
    % if the initial guess is all zeros, use the stiffness
    % matrix as the initial Jacobian to improve the guess
    u = [K+A+Q H'; H zeros(size(H,1))] \ [F+G; R];
else
    u=[u0 ;H'\((F+G)-(K+A+Q)*u0)];
end

if all(all(isfinite(u)))~=1
    error(message('pde:pdenonlin:InvalidInitGuess'))
end

global pdenonn pdenonTempPde pdenonp pdenont pdenoncoef
pdenonTempPde=self; pdenonp=p; pdenont=t; pdenoncoef=coefstruct;
pdenonn=np*self.PDESystemSize;
global pdenonK

r=compResidual(0,u);
K=pdenonK; % pdenonK may change after call to compResidual or numjac

if ischar(self.SolverOptions.ResidualNorm) % energy norm
    if strcmp(self.SolverOptions.ResidualNorm, 'energy')
        nr2=r'*K*r;
    else
        nr2=norm(r,self.SolverOptions.ResidualNorm)^2; % Frobenius norm - NOT DOCUMENTED!!!
    end
else
    nr2=(norm(r,self.SolverOptions.ResidualNorm)/length(r)^(1/self.SolverOptions.ResidualNorm))^2;
end
res=sqrt(nr2);
if strcmp(self.SolverOptions.ReportStatistics,'on')
    disp('Iteration     Residual     Step size  Jacobian: Full')
    fprintf('%4i%20.4e\n',0,res);
end

if (res > self.SolverOptions.ResidualTolerance)
    pseudocoefstruct = coefstruct;
    for i = 1:numel(pseudocoefstruct.Coefficients.c)
    pseudocoefstruct.Coefficients.c{i} = 1;
    pseudocoefstruct.Coefficients.a{i} = ones(N*N,1);
    pseudocoefstruct.Coefficients.f{i} = zeros(N,1);
    end
        
    [K, A, ~, ~, ~, H, ~] = assemblePdeMatrices(self, ndim, p, t, pseudocoefstruct, ones(size(u(1:pdenonn))));
    %                                         c, a         ,          f,        u
    % [a1,a2,a3,a4,a5,a6,a7]=assempde(b,p,e,t,0,ones(N*N,1), zeros(N,1), ones(np*N,1)); ???
    
    if ~any(any(H))
        H=sparse(0,size(H,2));
    end
    S=([K+A H';H zeros(size(H,1))]~=0);
    thresh = calcThresh(u, np*(self.PDESystemSize));
    pdeResidual = str2func('compResidual');
    [dfdu0,fac,G]=numjac(pdeResidual,0,u,r,thresh,[],0,S,[]);
    first_time=1;
end


iter=0;
res = zeros(self.SolverOptions.MaxIterations,1);
while (sqrt(nr2)> self.SolverOptions.ResidualTolerance)  % Gauss-Newton iteration loop
    iter=iter+1;
    if iter>self.SolverOptions.MaxIterations
        error(message('pde:pdenonlin:TooManyIter')) % *** TODO Change message ID
    end
    %J=K;
  
    if ~first_time
        thresh = calcThresh(u, np*N);
        [J,fac,G]=numjac(pdeResidual,0,u,r,thresh,fac,0,S,G);
    else
        J=dfdu0;
    end
    first_time=0;
    
    alfa=1;
    lastwarn('','');
    warnState = warning('off', 'MATLAB:singularMatrix');
    cor=-J\r;
    [~,wid] = lastwarn;
    if(strcmp(wid,'MATLAB:singularMatrix'))
        warning(warnState);
        error(message('pde:pdenonlin:SingularJacobian'))
    end
    warning(warnState);
    
    while(true)
        if alfa<self.SolverOptions.MinStep, error(message('pde:pdenonlin:SmallStepSize')), end
        tmp=u+alfa*cor;
        rr=compResidual(0,tmp);
        
        if ischar(self.SolverOptions.ResidualNorm) % energy norm
            if strcmp(self.SolverOptions.ResidualNorm, 'energy')
                nrr=rr'*K*rr;
            else
                nrr=norm(rr,self.SolverOptions.ResidualNorm)^2; % Frobenius norm - NOT DOCUMENTED!!!
            end
        else
            nrr=(norm(rr,self.SolverOptions.ResidualNorm)/length(rr)^(1/self.SolverOptions.ResidualNorm))^2;
        end
        
        if nr2-nrr < alfa/2*nr2
            alfa=alfa/2;
        else
            u=tmp;
            r=rr;
            nr2=nrr;
            K=pdenonK;
            if strcmp(self.SolverOptions.ReportStatistics,'on')
                fprintf('%4i%20.4e%12.7f\n',iter,sqrt(nr2),alfa);
            end
            break
        end
    end % while(true)
    res(iter)=sqrt(nr2);
end % END  - Gauss-Newton iteration loop

u=full(u(1:np*(self.PDESystemSize)));

clear global pdenonTempPde pdenonp  pdenont pdenoncoef
clear global pdenonK pdenonF pdenonn 
end

function [K, M, F, Q, G, H, R] = assemblePdeMatrices(thePde, ndims, p, t, coefstruct, u)

tempPde = pde.PDEModel(thePde.PDESystemSize);

if(ndims==2)
    [K, F] = formGlobalKF2D(tempPde, p, t, coefstruct,u,[]);
    M = formGlobalM2D(tempPde, p, t, coefstruct,u,[],'a');
elseif(ndims==3)
    [K, F] = formGlobalKF3D(tempPde, p, t, coefstruct,u,[]);
    M = formGlobalM3D(tempPde, p, t, coefstruct,u,[],'a');
end
% Impose BCs
[Q,G,H,R]=thePde.assembleBoundary(u);

end



function r=compResidual(~,u)
%compResidual Residual for nonlinear solver

%       Copyright 2015 The MathWorks, Inc.


global pdenonn pdenonTempPde pdenonp pdenont pdenoncoef
global pdenonK pdenonF

u=full(u);
ndims = size(pdenonp,1);
[K, M, F, Q, G, H, R] = assemblePdeMatrices(pdenonTempPde, ndims, pdenonp, pdenont, pdenoncoef, u(1:pdenonn));

pdenonK=K+M+Q;
pdenonF=F+G;
if any(any(H))
    if size(u,1)<(size(pdenonK,2)+size(H,1))
        u=[u;H'\(pdenonF-pdenonK*u)];
    end
    pdenonK=[pdenonK H';H zeros(size(H,1))];
    pdenonF=[pdenonF;R];
end
r=pdenonK*u-pdenonF;

end


function thresh = calcThresh(u, numDofs)
nu = size(u,1);
thresh=zeros(nu,1);
threshFact = 1e-6;
if all(u(1:numDofs)==0)
    thresh(1:numDofs)=sqrt(eps)*ones(numDofs,1);
else
    thresh(1:numDofs)=threshFact*max(abs(u(1:numDofs)))*ones(numDofs,1);
end
if numDofs<nu
    if all(u(numDofs+1:nu)==0)
        thresh(numDofs+1:nu)=sqrt(eps)*ones(nu-numDofs,1);
    else
        thresh(numDofs+1:nu)=threshFact*max(abs(u(numDofs+1:nu)))*ones(nu-numDofs,1);
    end
end
end