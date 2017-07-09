% Simulates a singlecell model for an entire bacterium
classdef Sim < handle

properties
  tstart;
  y;
end

methods
  function this= Sim(y0,t0)
    this.y= y0;
    if nargin==1, this.tstart= 0;
    else, this.tstart= t0; end
  end
  function extraAHL= step(this,AHL,tstop)
    this.y(6)= AHL;
    this.y= singlecell.run(this.y,[this.tstart tstop]);
    this.y= this.y(end,:);
    this.tstart= this.tstart + tstop;
    extraAHL= this.y(6)-AHL;
  end
end
end
