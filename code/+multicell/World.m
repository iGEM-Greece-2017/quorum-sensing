classdef World < handle
properties
  % Params: (can change on the fly)
  % - replicationPeriod: timesteps between bacterial replications
  % - singleBactModel:
  %   - 
  % - diffusionConstant: controls AHL diffusion strength
  params;
end
properties (SetAccess= private)
  timestep;
  worldSize;      % 2-dim
  hasBacterium;   % [Nx]x[Ny] {sparse<->full,logical}
  bacteriumState; % struct of many ([Nx]x[Ny] {sparse<->full})
  AHL;            % [Nx]x[Ny]
end

methods
  function this= World(worldSize, params)
    this.params= params;
    this.timestep= 0;
    this.worldSize= worldSize;
    this.hasBacterium= logical(sparse(worldSize,worldSize));
    this.AHL= sparse(worldSize,worldSize);
    % TODO: bacteriumState init
    % TODO: params (validate,fill)
  end
  
  function step(this)
  % Applies the 3 world transition functions: singleBactModel, AHL_diffuse, bact_lifecycle
    this.timestep= this.timestep+1;
    [this.bacteriumState, this.AHL]= multicell.singleBactModel(this.hasBacterium, ...
                          this.bacteriumState, this.AHL, this.params.singleBactModel);
    this.AHL= multicell.diffuse(this.AHL, this.params.diffusionConstant);
    % Don't replicate on every timestep, but each "replicationPeriod" timesteps
    if ~mod(this.timestep, this.params.replicationPeriod)
      this.hasBacterium= multicell.lifecycle(this.hasBacterium);
    end
    
    % Transform state representation between sparse and full based on bacterial coverage
    bacterialCoverage= nnz(this.hasBacterium)/numel(this.hasBacterium);
    if bacterialCoverage > this.params.sparsity.high && issparse(this.hasBacterium)
      this.hasBacterium= full(this.hasBacterium);
      % TODO: this.bacteriumState
    elseif bacterialCoverage < this.params.sparsity.low && ~issparse(this.hasBacterium)
      this.hasBacterium= sparse(this.hasBacterium);
      % TODO: this.bacteriumState
    end
  end
  
  function [timestep,worldSize,hasBacterium,AHL,bacteriumState]= getState(this)
    timestep= this.timestep;
    worldSize= this.worldSize;
    hasBacterium= this.hasBacterium;
    AHL= this.AHL;
    bacteriumState= this.bacteriumState;
  end
end
end
