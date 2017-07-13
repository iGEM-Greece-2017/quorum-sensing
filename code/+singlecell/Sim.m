% Simulates a singlecell model for an entire bacterium
classdef Sim < handle

properties
  stateHistory;
end

methods
  function this= Sim(y0,t0)
    this.stateHistory.y(:,1)= y0;
    this.stateHistory.t(1)= 0;
    this.stateHistory.latest= 1;
    if nargin==2, this.stateHistory.t(1)= t0; end
  end
  function ahlProd= step(this,AHL,tstop)
    ahlProd= 0;
    latest= this.rollback(tstop);
    if latest>0
      y= this.stateHistory.y(:,latest);
      tstart= this.stateHistory.t(latest);
      y= singlecell.run(y,[tstart tstop])';
      ahlProd= (y(6,end) - AHL) ./ (tstop - tstart);
      this.stateHistory.y(:,latest+1)= y(:,end);
      this.stateHistory.y(6,latest+1)= AHL;     % Enforce pde AHL value
      this.stateHistory.t(latest+1)= tstop;
      this.stateHistory.latest= latest+1;
    end
  end
end

methods (Access= private)
  function latest= rollback(this,tstop)
    latest= 0;
    for i= this.stateHistory.latest:-1:1
      if this.stateHistory.t(i) < tstop
        latest= i;
        return;
      end
    end
  end
end
end
