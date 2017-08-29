function f=firstOrderODEf(t,u)
%firstOrderODEf - Residue function for first-order ODE system
%Computes residue of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015 The MathWorks, Inc.

global femodel
%{ @@@ Custom @@@ %%
global enableSinglecellEq;
global bactNodes;
global bactNodeCoeffs;
nBact= length(bactNodes);

femodelLoc= femodel;
enableSinglecellEqLoc= enableSinglecellEq;
bactNodesLoc= bactNodes;
bactNodeCoeffsLoc= bactNodeCoeffs;
%} @@@ Custom @@@ %%

f= zeros(size(u));
newK= cell(size(u,2),1); newF= cell(size(u,2),1); A= cell(size(u,2),1);
parfor i= 1:size(u,2)
  returnFlag= false;
  
  if ~(femodelLoc.vq || femodelLoc.vg || femodelLoc.vh || femodelLoc.vr || ...
          femodelLoc.vc || femodelLoc.va || femodelLoc.vf)
      if size(u,2)>1, error('This codepath can''t be used with vectorized=''on'''); end
      if enableSinglecellEqLoc
        f(:,i)=[-femodelLoc.K*u(1:end-nBact*8)+femodelLoc.F; zeros(nBact*8,1)];
      else, f(:,i)=-femodelLoc.K*u+femodelLoc.F;
      end
      returnFlag= true;
  end
  
  if ~returnFlag
    % make uFull
    if(femodelLoc.numConstrainedEqns == femodelLoc.totalNumEqns)
        if enableSinglecellEqLoc, uFull= u(1:end-nBact*8,i);
        else, uFull= u(:,i); end
    else
        if size(u,2)>1, error('This codepath can''t be used with vectorized=''on'''); end
        if enableSinglecellEqLoc, uFull= femodelLoc.B*u(1:end-nBact*8) + femodelLoc.ud;
        else, uFull= femodelLoc.B*u + femodelLoc.ud; end
    end
    
    % [K,F,A]= formGlobalKF2D(uFull)
    if femodelLoc.vc || femodelLoc.va || femodelLoc.vf
      if(femodelLoc.nrp==2)
          [newK{i}, newF{i}]= formGlobalKF2D(femodelLoc.thePde, femodelLoc.p, femodelLoc.t, femodelLoc.coefstruct,uFull,t);
          A{i}= formGlobalM2D(femodelLoc.thePde, femodelLoc.p, femodelLoc.t, femodelLoc.coefstruct,uFull,t,'a');
      elseif(femodelLoc.nrp==3)
          error('Not implemented!');
      end
    end

    if ~(femodelLoc.vh || femodelLoc.vr)
        % neither H nor R are functions of t or u
        K = newK{i} + femodelLoc.A + femodelLoc.Q;
        F = femodelLoc.B'*(newF{i} + femodelLoc.G - K*femodelLoc.ud);
        K = femodelLoc.B'*K*femodelLoc.B;
    else
        error('Not implemented!');
    end

  %% @@@ Custom @@@ %%
    if ~enableSinglecellEqLoc
      f(:,i)= -K*u(:,i)+F;
      returnFlag= true;
    end
    if ~returnFlag
      % Initial calculation for d(n_i)/dt
      fi= [-K*u(1:end-nBact*8,i)+F; zeros(nBact*8,1)];

      for b= 1:nBact
        ahl= mean(u(bactNodesLoc{b},i));
        yIdx= (nBact-b+1)*8-1;    % where the y variables for this bacterium start (relative to end)
        yBact= [u(end-yIdx: end-yIdx+4,i); ahl; u(end-yIdx+5: end-yIdx+7,i)];   % y variables for this bacterium

        dy= singlecell.model(1,yBact);
        ahlProd= dy(6).*bactNodeCoeffsLoc{b};

        % Note:
        %   The cylindrical symmetry allows for many bacteria that sit on a cylindrical ring
        % of set radius to be simulated by a single set of singlecell equations. Their total
        % effect is incorporated by multiplying their number with their output (dy_6)
        fi(bactNodesLoc{b})= fi(bactNodesLoc{b}) + ahlProd;

        fi(end-yIdx: end-yIdx+7)= dy([1:5,7:9]);
      end

      % Don't allow nodes to stay negative
      %fi(u<=0 & f<=0)= -fi(u<=0 & f<=0);%-u(u<=0 & f<=0);
      f(:,i)= fi;
    end
  end
end

femodelLoc.A= A{end}; % TEMP - Change the femodelLoc to have a A matrix later.
femodelLoc.F= newF{end};
femodelLoc.K= newK{end};

femodel= femodelLoc;
bactNodes= bactNodesLoc;
bactNodeCoeffs= bactNodeCoeffsLoc;
