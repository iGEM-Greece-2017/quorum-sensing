% Holds a vector of singlecell models for the bacteria that are going to be simulated and, 
% given a point's (x,y) coordinates, routes the AHL production queries to the corresponding model
classdef BacteriaTable < handle

properties
  cellModels;   % Vector of singlecell models
end

methods
  function this= BacteriaTable(nBact)
    y0= [1.6E-7;0;0;0;0;0;0;0;0];
    t0= 0;
    % TODO: initialize singlecell.Sim correctly (check all objects)
    cellModels(nBact)= singlecell.Sim(y0,t0);
    this.cellModels= cellModels;
  end
  % This function is meant to be called by the pde solver on all nodes
  % inside a bacterium
  function extraAHL= AHLProd(this,region,state)
    n= size(region.x);    % How many points
    if state.time==0, extraAHL= zeros(n); return; end
    cell= this.routeCoordCell(region);
    %region
    state
    region.subdomain
    %extraAHL= this.cellModels(cell).step(state.u,state.time);
    extraAHL= 1e-6;
  end
end

methods (Access= private)
  % Maps coordinates to singlecell models and creates new models, if required
  function cellID= routeCoordCell(this,coord)
    % TODO
    cellID= 1;
  end
end
end
