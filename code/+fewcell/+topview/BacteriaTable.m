% Holds a vector of singlecell models for the bacteria that are going to be simulated and, 
% given a point's (x,y) coordinates, routes the AHL production queries to the corresponding model
classdef BacteriaTable < handle

properties
  cellModels;   % Vector of singlecell models
end

methods
  function this= BacteriaTable(nBact,y0,t0)
    % TODO: initialize singlecell.Sim correctly (check all objects)
    for b= 1:nBact+1, cellModels(b)= singlecell.Sim(y0,t0); end
    this.cellModels= cellModels;
  end
  
  % This function is meant to be called by the pde solver on all nodes
  % inside a bacterium
  function extraAHL= AHLProd(this,region,state)
    if any(isnan(state.u)), extraAHL= nan(size(state.u)); return; end
    if state.time < 1e-6, extraAHL= zeros(size(state.u)); return; end
    
    cellIdx= this.routeCoordCell(region);
    extraAHL= zeros(size(region.x));
    extraAHL(cellIdx==0)= 0;
    for i= unique(cellIdx)
      if i>0
        % Pass state.time to singlecell and let it match the given time with the nearest time
        % asked in the past, then load that state, erase all future times and roll forward
        % to the given time
        extraAHL(cellIdx==i)= this.cellModels(i).step(...
          mean(state.u(cellIdx==i)),state.time) ./ sum(cellIdx==i);
      end
    end
  end
end

methods (Access= private)
  % Maps coordinates to singlecell models and creates new models, if required
  function cellIdx= routeCoordCell(this,region)
    % Assigns each element to a cell.
    % Output:
    % - idx {size(region.x)}: cell to which each element is assigned (0 => no cell)
    
    % TODO: Only for 1 bact!
    cellIdx= abs(round(region.subdomain*1e300*1e23)-1);
  end
end
end
