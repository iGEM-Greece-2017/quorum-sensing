classdef World < handle
properties
  % Params: (can change on the fly)
  % - replicationPeriod: timesteps between bacterial replications
  % - singleBactModel:
  %   - 
  % - diffusion: diffusion constants
  %   - D_AHL
  %   - C_agar
  %   - t_s {derived}
  params;
end
properties %(SetAccess= private)
  timestep;
  worldSize;      % 2-dim
  hasBacterium;   % [Nx]x[Ny] %{sparse<->full,logical}
  bacteriumState; % [Nx]x[Ny] cellarray of [9] elements
  AHL;            % [Nx]x[Ny]
  producedAHL;    % [Nx]x[Ny]: buffer for the AHL produced in batch from the singlebact model and released gradually
  hasGPU;
end

methods
  function this= World(worldSize, randPrc, params)
    this.hasGPU= gpuDeviceCount>0;
    this.params= multicell.World.fillParams(params);
    this.timestep= 0;
    this.worldSize= worldSize;
    this.AHL= zeros(worldSize);
    this.producedAHL= zeros(worldSize);
    if this.hasGPU, this.AHL= gpuArray(this.AHL); end
    % Initial bacterial populations: only in the middle half of the board
    this.hasBacterium= false(worldSize);
    bactPlacementSize= size(this.hasBacterium(ceil(worldSize(1)/4)+1:worldSize(1)-floor(worldSize(1)/4), ...
      ceil(worldSize(2)/4)+1:worldSize(2)-floor(worldSize(2)/4)));
    this.hasBacterium(ceil(end/4)+1:end-floor(end/4), ceil(end/4)+1:end-floor(end/4))= ...
      rand(bactPlacementSize) < prod(worldSize)/prod(bactPlacementSize)*randPrc;
    this.bacteriumState= cellfun(@(x) [this.params.initDNA; zeros(8,1)], ...
      cell(worldSize), 'uniformoutput',0);
  end
  
  function step(this)
  % Applies the 3 world transition functions: singleBactModel, AHL_diffuse, bact_lifecycle
    this.AHL= multicell.diffuse(this.AHL,this.params.dt,this.params.r, this.params.diffusion);
    % Adds one part of the AHL produced in the prev singleBact model
    this.AHL(this.hasBacterium)= this.AHL(this.hasBacterium) + this.producedAHL(this.hasBacterium);
    % Run single bacterium model only once every so many iterations
    if ~mod(this.timestep, this.params.period.singleBactModel)
      if this.hasGPU, this.AHL= gather(this.AHL); end
      bactTimestep= floor(this.timestep/this.params.period.singleBactModel);
      [this.bacteriumState, this.producedAHL]= multicell.singleBactModel(this.hasBacterium, ...
                    this.bacteriumState, this.AHL, bactTimestep, this.params.singleBactModel);
      this.producedAHL= this.producedAHL / this.params.period.singleBactModel;
      if this.hasGPU, this.AHL= gpuArray(this.AHL); end
    end
    % Don't replicate on every timestep, but once each "replicationPeriod" timesteps
    if ~mod(this.timestep, this.params.period.replication)
      this.hasBacterium= multicell.lifecycle(this.hasBacterium);
    
      %{
      % Transform state representation between sparse and full based on bacterial coverage
      bacterialCoverage= nnz(this.hasBacterium)/numel(this.hasBacterium);
      if bacterialCoverage > this.params.sparsity.high && issparse(this.hasBacterium)
        this.hasBacterium= full(this.hasBacterium);
        % TODO: this.bacteriumState
      elseif bacterialCoverage < this.params.sparsity.low && ~issparse(this.hasBacterium)
        this.hasBacterium= sparse(this.hasBacterium);
        % TODO: this.bacteriumState
      end
      %}
    end
    
    this.timestep= this.timestep+1;
    %% Debug
    % Total bacteria, total AHL
    %fprintf('[World.step]: total bacteria: %d \n',sum(sum(this.hasBacterium)));
    %fprintf('[World.step]: total AHL: %g \n', sum(sum(this.AHL)));
  end
  
  function [timestep,worldSize,hasBacterium,AHL,bacteriumState]= getState(this)
    timestep= this.timestep;
    worldSize= this.worldSize;
    hasBacterium= this.hasBacterium;
    AHL= this.AHL;
    bacteriumState= this.bacteriumState;
  end
end
methods (Static)
  function params= fillParams(params)
    % Existence assertions
    assert(isfield(params,'dt'));
    assert(isfield(params,'r'));
    assert(isfield(params,'initDNA'));
    
    assert(isfield(params,'period'));
    assert(isfield(params.period,'singleBactModel'));
    assert(isfield(params.period,'replication'));
    
    assert(isfield(params,'diffusion'));
    assert(isfield(params.diffusion,'D_AHL'));
    assert(isfield(params.diffusion,'C_agar'));
    
    % Fill in the blanks
    params.diffusion.t_s= (2*params.r)^2 / (4*params.diffusion.C_agar*params.diffusion.D_AHL);
    params.singleBactModel.dt= params.dt * params.period.singleBactModel;
  end
end
end
